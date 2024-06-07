#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good BAM files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/bam/good
DIR_BASIC="$DIR_OUT/basic.bam"

# ------------------------------------------------------------------------------
# SAM files
# ------------------------------------------------------------------------------

log "Creating valid SAM file"
samtools view -h "$DIR_BASIC" > "$DIR_OUT/basic.sam"
validate "ok"

log "Creating compressed SAM files (bgzip)"
cp "$DIR_OUT/basic.sam" "$DIR_OUT/compressed.sam"
bgzip -f "$DIR_OUT/compressed.sam"
validate "ok"

log "Creating indexed SAM files (CSI, TBI)"
cp "$DIR_OUT/compressed.sam.gz" "$DIR_OUT/indexed_csi.sam.gz"
cp "$DIR_OUT/compressed.sam.gz" "$DIR_OUT/indexed_tbi.sam.gz"
tabix --csi -p sam "$DIR_OUT/indexed_csi.sam.gz"
tabix -p sam "$DIR_OUT/indexed_tbi.sam.gz"
validate "ok"

# ------------------------------------------------------------------------------
# Indexed
# ------------------------------------------------------------------------------

log "Creating indexed BAM files (BAI, CSI)"
cp "$DIR_BASIC" "$DIR_OUT/indexed_bai.bam"
cp "$DIR_BASIC" "$DIR_OUT/indexed_csi.bam"
samtools index --bai "$DIR_OUT/indexed_bai.bam"
samtools index --csi "$DIR_OUT/indexed_csi.bam"
validate "ok"

# ------------------------------------------------------------------------------
# No mapping information
# ------------------------------------------------------------------------------

log "Creating BAM with no mapped reads"
samtools reset "$DIR_BASIC" -O BAM > "$DIR_OUT/no_mapped_reads.bam"
validate "$(samtools quickcheck "$DIR_OUT/no_mapped_reads.bam" 2>&1 | grep "no targets in header")"


# ==============================================================================
# Bad BAM files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/bam/bad
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Truncated
# ------------------------------------------------------------------------------

log "Creating truncated BAM file"
head -c 10000 "$DIR_BASIC" > "$DIR_OUT/truncated.bam"
validate "$(samtools view "$DIR_OUT/truncated.bam" 2>&1 | grep "EOF marker is absent. The input is probably truncated")"

# ------------------------------------------------------------------------------
# BAI index is older than the BAM file (might needs reindexing to reflect latest BAM file contents)
# ------------------------------------------------------------------------------

log "Creating BAI index that is older than the BAM file"
cp "$DIR_BASIC" "$DIR_OUT/bai_older_than_data.bam"
samtools index "$DIR_OUT/bai_older_than_data.bam"
sleep 1  # wait before run touch to make sure the timestamp is different enough
touch "$DIR_OUT/bai_older_than_data.bam"
validate "$(samtools idxstats "$DIR_OUT/bai_older_than_data.bam" 2>&1 | grep "The index file is older than the data file")"
