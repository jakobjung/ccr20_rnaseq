#!/bin/bash
# Trimming, mapping and counting pipeline for the CCR20 PNA RNA-seq project.
# Paired-end reads, single reference (E. coli 536 / CP000247.1, a UPEC strain
# carrying the colibactin (pks) island). Designed to run on the Sanger farm
# (LSF/bsub) via one job per sample.
#
# Expected layout (relative to this script, i.e. this repo checked out on
# /lustre/.../rnaseq_upec/):
#   ../data/fastq/X208SC26074449-Z01-F001/01.RawData/<SAMPLE>/<SAMPLE>_*_L*_{1,2}.fq.gz
#   ../data/reference_sequences/ecoli536.fasta
#   ../data/reference_sequences/ecoli536_sRNAs_modified.gff3
#   ../data/reference_sequences/adapters.fa
#
# Each sample was sequenced on two flowcells/lanes (top-up run), so the raw
# R1/R2 files for a sample are concatenated before trimming.
#
# Usage:
#   bash trimm_map_BB.sh submit          # bsub one trim+map job per sample
#   bash trimm_map_BB.sh run <SAMPLE>    # trim+map a single sample (called by bsub)
#   bash trimm_map_BB.sh counts          # bsub featureCounts once all samples are mapped
#
# After trim+map jobs finish, check bhist/bjobs, then run "counts".

set -euo pipefail

# ---- paths -----------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="$SCRIPT_DIR/../data"
RAWDATA="$PROJECT/fastq/X208SC26074449-Z01-F001/01.RawData"
REF_FASTA="$PROJECT/reference_sequences/ecoli536.fasta"
REF_GFF="$PROJECT/reference_sequences/ecoli536_sRNAs_modified.gff3"
ADAPTERS="$PROJECT/reference_sequences/adapters.fa"
LIBDIR="$PROJECT/libs"
ALIGNDIR="$PROJECT/rna_align"
LOGDIR="$SCRIPT_DIR/logs"

# ---- farm / LSF settings (edit to match what's available on farm22) -------
QUEUE="normal"
CORES=4
MEM_MB=8000            # per bsub job, in MB
JAVA_XMX="6g"           # bbduk/bbmap JVM heap, keep below MEM_MB. Must be set
                         # explicitly: bbmap.sh otherwise auto-sizes its heap
                         # from the *node's* total RAM (up to ~2TB on some
                         # farm22 hosts), not the LSF job's memory reservation,
                         # and gets silently killed for exceeding -M.
GROUP="team377f"        # LSF fairshare group (from `show_my_lsf_groups`); leave empty to omit -G

BSUB_GROUP_OPTS=()
[ -n "$GROUP" ] && BSUB_GROUP_OPTS=(-G "$GROUP")

# Load these before running (module names as available on farm22):
#   module load bbtools/39.01 samtools/1.21 subread/2.0.6--he4a0461_0
# NB subread must be >=2.0.3 for the --countReadPairs flag used in run_counts()
# below — do NOT use subread/1.4.5-p1, it predates that flag and featureCounts
# will error out on it.

list_samples(){
    find "$RAWDATA" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
}

submit_all(){
    mkdir -p "$LOGDIR"
    for SAMPLE in $(list_samples); do
        bsub -q "$QUEUE" -n "$CORES" "${BSUB_GROUP_OPTS[@]}" \
             -R "span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" -M "$MEM_MB" \
             -J "map_${SAMPLE}" \
             -o "$LOGDIR/${SAMPLE}.out" -e "$LOGDIR/${SAMPLE}.err" \
             "bash '$SCRIPT_DIR/trimm_map_BB.sh' run '$SAMPLE'"
    done
}

run_sample(){
    SAMPLE="$1"
    SAMPLEDIR="$RAWDATA/$SAMPLE"
    if [ ! -d "$SAMPLEDIR" ]; then
        echo "No such sample directory: $SAMPLEDIR" >&2
        exit 1
    fi

    mkdir -p "$LIBDIR" "$ALIGNDIR"

    R1_CAT="$LIBDIR/${SAMPLE}_1.fq.gz"
    R2_CAT="$LIBDIR/${SAMPLE}_2.fq.gz"
    R1_TRIM="$LIBDIR/${SAMPLE}_1_trimmed.fq.gz"
    R2_TRIM="$LIBDIR/${SAMPLE}_2_trimmed.fq.gz"

    echo "[$SAMPLE] concatenating reads across flowcells/lanes"
    cat "$SAMPLEDIR"/*_1.fq.gz > "$R1_CAT"
    cat "$SAMPLEDIR"/*_2.fq.gz > "$R2_CAT"

    echo "[$SAMPLE] trimming adapters/quality with bbduk"
    # NB this is a first quick-and-dirty pass with generic settings; check the
    # fastqc output and adjust ktrim/qtrim/trimq if adapter or quality content
    # looks unusual for this library prep.
    bbduk.sh in1="$R1_CAT" in2="$R2_CAT" \
             ref="$ADAPTERS" -Xmx${JAVA_XMX} t="$CORES" \
             out1="$R1_TRIM" out2="$R2_TRIM" \
             ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10

    echo "[$SAMPLE] mapping to ecoli536 with bbmap"
    # nodisk=t: build the index in memory instead of caching it on disk under
    # ./ref/ (relative to CWD) — with all samples run as concurrent LSF jobs
    # sharing the same CWD, an on-disk cache would race across jobs. The
    # genome is tiny (~5Mb) so rebuilding per job costs a couple of seconds.
    bbmap.sh in1="$R1_TRIM" in2="$R2_TRIM" \
             ref="$REF_FASTA" t="$CORES" trimreaddescription=t k=12 -Xmx${JAVA_XMX} nodisk=t \
             outm="$ALIGNDIR/${SAMPLE}.sam"

    echo "[$SAMPLE] sorting & indexing bam"
    samtools sort -O BAM -@ "$CORES" "$ALIGNDIR/${SAMPLE}.sam" > "$ALIGNDIR/${SAMPLE}.bam"
    samtools index "$ALIGNDIR/${SAMPLE}.bam"
    rm "$ALIGNDIR/${SAMPLE}.sam"

    echo "[$SAMPLE] fastqc on trimmed reads"
    fastqc -q "$R1_TRIM" "$R2_TRIM" -o "$LIBDIR" || echo "[$SAMPLE] fastqc failed, continuing"

    echo "[$SAMPLE] done"
}

build_saf(){
    # featureCounts doesn't reliably parse GFF3's key=value attribute syntax,
    # so build a plain SAF file (GeneID/Chr/Start/End/Strand) straight from
    # the CDS/sRNA/tRNA/rRNA rows' locus_tag= field instead of feeding it the
    # gff3 directly.
    SAF="$ALIGNDIR/ecoli536_cds_sRNA_tRNA_rRNA.saf"
    mkdir -p "$ALIGNDIR"
    {
        printf 'GeneID\tChr\tStart\tEnd\tStrand\n'
        awk -F'\t' '$3=="CDS" || $3=="sRNA" || $3=="tRNA" || $3=="rRNA" {
            match($9, /locus_tag=[^;]+/)
            tag = substr($9, RSTART+10, RLENGTH-10)
            print tag "\t" $1 "\t" $4 "\t" $5 "\t" $7
        }' "$REF_GFF"
    } > "$SAF"
    echo "$SAF"
}

run_counts(){
    mkdir -p "$LOGDIR"
    SAF="$(build_saf)"
    bsub -q "$QUEUE" -n "$CORES" "${BSUB_GROUP_OPTS[@]}" \
         -R "span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" -M "$MEM_MB" \
         -J featurecounts -o "$LOGDIR/featurecounts.out" -e "$LOGDIR/featurecounts.err" \
         "featureCounts -T $CORES -p --countReadPairs -F SAF \
             -a '$SAF' -o '$ALIGNDIR/counttable_ccr20.txt' '$ALIGNDIR'/*.bam"
}

case "${1:-}" in
    submit) submit_all ;;
    run)    run_sample "${2:?sample name required}" ;;
    counts) run_counts ;;
    *)
        echo "Usage: $0 {submit|run <SAMPLE>|counts}"
        exit 1
        ;;
esac
