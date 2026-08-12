#!/bin/bash
# Trimming, mapping and counting pipeline for the CCR20 PNA RNA-seq project.
# Paired-end reads, mapped against a single reference (E. coli 536 /
# CP000247.1, a UPEC strain carrying the colibactin (pks) island).
#
# Run the whole thing as a single job, e.g.:
#   bsub -q long -n 4 -M 10000 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" \
#        -o scripts/stdout.%J -e scripts/stderr.%J bash scripts/trimm_map_BB.sh
# (normal's 12h run limit may be tight for 42 samples run sequentially; long
# gives 48h. Add -G <group> if your team uses LSF fairshare groups.)
#
# Load these first (module names as available on farm22):
#   module load bbtools/39.01 samtools/1.21 subread/2.0.6--he4a0461_0
# NB subread must be >=2.0.3 for the --countReadPairs flag used below — do not
# use subread/1.4.5-p1 or 2.0.1, they predate that flag.

main(){
    # resolve relative to this script's own location, not the caller's CWD,
    # so it works whether you submit it from scripts/ or the repo root.
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    PROJECT=$SCRIPT_DIR/../data
    echo "Start trimming"
    rename_trim_rna_libs
    echo "Trimming done. Start mapping"
    align_rna_reads_genome
    echo "Finished mapping. Start counting"
    count_features
}

rename_trim_rna_libs(){
    mkdir -pv $PROJECT/libs
    RAWDATA=$PROJECT/fastq/X208SC26074449-Z01-F001/01.RawData
    for SAMPLEDIR in $RAWDATA/*/
    do
        SAMPLE=$(basename $SAMPLEDIR)
        echo "$SAMPLE starts trimming nowwwwwww"
        # each sample was sequenced across two flowcells/lanes, concatenate first:
        cat $SAMPLEDIR/*_1.fq.gz > $PROJECT/libs/${SAMPLE}_1.fq.gz
        cat $SAMPLEDIR/*_2.fq.gz > $PROJECT/libs/${SAMPLE}_2.fq.gz
        # bbduk trims low quality bases and removes adapters:
        bbduk.sh in1=$PROJECT/libs/${SAMPLE}_1.fq.gz in2=$PROJECT/libs/${SAMPLE}_2.fq.gz \
                 ref=$PROJECT/reference_sequences/adapters.fa -Xmx6g t=4 \
                 out1=$PROJECT/libs/${SAMPLE}_1_trimmed.fq.gz out2=$PROJECT/libs/${SAMPLE}_2_trimmed.fq.gz \
                 ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10
    done
}

align_rna_reads_genome(){
    mkdir -p $PROJECT/rna_align
    DIR=$PROJECT/rna_align
    for i in $(ls $PROJECT/libs/*_1_trimmed.fq.gz)
    do
        NAME=${i##*/}
        NAME=${NAME%_1_trimmed.fq.gz}
        echo "Starting mapping for sample: $NAME"
        # -Xmx must be set explicitly: bbmap.sh otherwise auto-sizes its heap
        # from the node's total RAM (not the LSF job's -M reservation) and
        # gets silently killed for exceeding it.
        bbmap.sh in1=$i in2=${i%_1_trimmed.fq.gz}_2_trimmed.fq.gz trimreaddescription=t t=4 \
                  ref=$PROJECT/reference_sequences/ecoli536.fasta k=13 -Xmx6g outm=$DIR/$NAME.sam
        # sort sam file, create BAM file:
        samtools sort -O BAM -@ 4 $DIR/$NAME.sam > $DIR/$NAME.bam
        # remove sam file: (not actually needed)
        rm $DIR/$NAME.sam
        # index BAM file:
        samtools index $DIR/$NAME.bam
        fastqc $i
    done
}

count_features(){
    # featureCounts doesn't reliably parse GFF3's key=value attribute syntax
    # (only GTF/GFF2's key "value";), so build a plain SAF file
    # (GeneID/Chr/Start/End/Strand) from the gff3's CDS/sRNA/tRNA/rRNA rows'
    # locus_tag= field instead of passing the gff3 to featureCounts directly.
    SAF=$PROJECT/reference_sequences/ecoli536_cds_sRNA_tRNA_rRNA.saf
    { echo -e "GeneID\tChr\tStart\tEnd\tStrand"
      awk -F'\t' '$3=="CDS" || $3=="sRNA" || $3=="tRNA" || $3=="rRNA" {
          match($9, /locus_tag=[^;]+/)
          tag = substr($9, RSTART+10, RLENGTH-10)
          print tag "\t" $1 "\t" $4 "\t" $5 "\t" $7
      }' $PROJECT/reference_sequences/ecoli536_sRNAs_modified.gff3
    } > $SAF

    featureCounts -T 4 -p --countReadPairs -F SAF \
                  -a $SAF -o $PROJECT/rna_align/counttable_ccr20.txt \
                  $PROJECT/rna_align/*.bam
}

main
