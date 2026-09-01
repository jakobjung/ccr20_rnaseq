#!/bin/bash
# Gathers the raw (untrimmed) paired FASTQ files for the 21 4h samples used in
# the analysis into one folder, ready for GEO/SRA upload, and writes their
# MD5 checksums (GEO requires these).
#
# Run this ON THE FARM (where data/fastq/.../01.RawData/ actually lives), from
# anywhere - paths resolve relative to this script's own location, same as
# trimm_map_BB.sh:
#   bash scripts/gather_raw_fastq_4h.sh
#
# For each sample it prefers data/libs/<label>_1.fq.gz / _2.fq.gz - the
# lane-merged-but-untrimmed reads trimm_map_BB.sh already produced as bbduk's
# input (still genuinely "raw" - only the lanes were concatenated, nothing was
# trimmed/filtered). If those were since cleaned up, it falls back to
# concatenating fresh from the original per-lane delivery in 01.RawData/.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT=$SCRIPT_DIR/../data
RAWDATA=$PROJECT/fastq/X208SC26074449-Z01-F001/01.RawData
SAMPLE_TABLE=$PROJECT/20260709_RNAseq_sample_table.csv
DEST=$SCRIPT_DIR/../GEO_submission/raw_fastq_4h

mkdir -p "$DEST"

# columns in the sample table: label,RNA_conc_ng_ul,PNA,strain,clone,treatment_time_h,unique_identifier
SAMPLES=$(awk -F',' 'NR>1 && $6==4 {print $1","$7}' "$SAMPLE_TABLE")
echo "4h samples to gather: $(echo "$SAMPLES" | wc -l)"

while IFS=',' read -r LABEL SAMPLE_ID; do
  DEST_1="$DEST/${SAMPLE_ID}_R1.fastq.gz"
  DEST_2="$DEST/${SAMPLE_ID}_R2.fastq.gz"

  if [ -f "$PROJECT/libs/${LABEL}_1.fq.gz" ] && [ -f "$PROJECT/libs/${LABEL}_2.fq.gz" ]; then
    echo "$SAMPLE_ID ($LABEL): copying from data/libs/ (already lane-merged, pre-trim)"
    cp "$PROJECT/libs/${LABEL}_1.fq.gz" "$DEST_1"
    cp "$PROJECT/libs/${LABEL}_2.fq.gz" "$DEST_2"
  elif [ -d "$RAWDATA/$LABEL" ]; then
    echo "$SAMPLE_ID ($LABEL): data/libs/ copy not found - concatenating lanes fresh from 01.RawData/"
    cat "$RAWDATA/$LABEL"/*_1.fq.gz > "$DEST_1"
    cat "$RAWDATA/$LABEL"/*_2.fq.gz > "$DEST_2"
  else
    echo "WARNING: no raw data found for sample '$LABEL' ($SAMPLE_ID) in data/libs/ or 01.RawData/ - skipping" >&2
    continue
  fi
done <<< "$SAMPLES"

echo "Computing MD5 checksums..."
(cd "$DEST" && md5sum *.fastq.gz > md5sums_4h.txt)

echo "Done. $(ls "$DEST"/*.fastq.gz | wc -l) files gathered in: $DEST"
echo "Checksums: $DEST/md5sums_4h.txt"
