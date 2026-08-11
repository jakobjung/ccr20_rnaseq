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
BBDUK_XMX="6g"          # bbduk/bbmap JVM heap, keep below MEM_MB

# TODO: load whatever gives you bbduk.sh/bbmap.sh, samtools and featureCounts
# on farm22, e.g.:
#   module load bbmap/38.84 samtools/1.12 subread/2.0.1
# or activate a conda env with these tools installed.

list_samples(){
    find "$RAWDATA" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
}

submit_all(){
    mkdir -p "$LOGDIR"
    for SAMPLE in $(list_samples); do
        bsub -q "$QUEUE" -n "$CORES" \
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
             ref="$ADAPTERS" -Xmx${BBDUK_XMX} t="$CORES" \
             out1="$R1_TRIM" out2="$R2_TRIM" \
             ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10

    echo "[$SAMPLE] mapping to ecoli536 with bbmap"
    bbmap.sh in1="$R1_TRIM" in2="$R2_TRIM" \
             ref="$REF_FASTA" t="$CORES" trimreaddescription=t k=12 \
             outm="$ALIGNDIR/${SAMPLE}.sam"

    echo "[$SAMPLE] sorting & indexing bam"
    samtools sort -O BAM -@ "$CORES" "$ALIGNDIR/${SAMPLE}.sam" > "$ALIGNDIR/${SAMPLE}.bam"
    samtools index "$ALIGNDIR/${SAMPLE}.bam"
    rm "$ALIGNDIR/${SAMPLE}.sam"

    echo "[$SAMPLE] fastqc on trimmed reads"
    fastqc -q "$R1_TRIM" "$R2_TRIM" -o "$LIBDIR" || echo "[$SAMPLE] fastqc failed, continuing"

    echo "[$SAMPLE] done"
}

run_counts(){
    mkdir -p "$LOGDIR"
    # locus_tag is used consistently for CDS/sRNA/tRNA/rRNA in this gff, so a
    # single featureCounts call covers all feature types at once.
    bsub -q "$QUEUE" -n "$CORES" \
         -R "span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" -M "$MEM_MB" \
         -J featurecounts -o "$LOGDIR/featurecounts.out" -e "$LOGDIR/featurecounts.err" \
         "featureCounts -T $CORES -p --countReadPairs -t CDS,sRNA,tRNA,rRNA -g locus_tag \
             -a '$REF_GFF' -o '$ALIGNDIR/counttable_ccr20.txt' '$ALIGNDIR'/*.bam"
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
