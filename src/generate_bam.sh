#!/bin/bash

. ./src/lib.sh || exit

# ==============================================================================
# Good BAM files
# ==============================================================================

URL="https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
DIR_OUT=data/bam/good

# ------------------------------------------------------------------------------
# Small BAM file
# ------------------------------------------------------------------------------

log "Creating valid BAM file"
samtools view -b -h "$URL" 11:82365011-82366010 > "$DIR_OUT/basic.bam"
validate $(samtools quickcheck -vvv "$DIR_OUT/basic.bam" 2>&1 | grep "has good EOF block")


# ==============================================================================
# Bad BAM files
# ==============================================================================

DIR_OUT=data/bam/bad

# ------------------------------------------------------------------------------
# Truncated BAM file
# ------------------------------------------------------------------------------

log "Creating truncated BAM file"
curl -s -H "Range: bytes=0-25000" "$URL" > "$DIR_OUT/truncated.bam"

# Samtools outputs a warning about an absent EOF marker
validate $(samtools view "$DIR_OUT/truncated.bam" 2>&1 | grep "EOF marker is absent. The input is probably truncated")

# ------------------------------------------------------------------------------
# BAI index is older than the BAM file (might needs reindexing to reflect latest BAM file contents)
# ------------------------------------------------------------------------------

log "Creating BAI index that is older than the BAM file"
samtools view -b -h "$URL" 11:82365011-82366010 > "$DIR_OUT/bai_older_than_data.bam"
samtools index "$DIR_OUT/bai_older_than_data.bam"
sleep 1  # wait before run touch to make sure the timestamp is different enough
touch "$DIR_OUT/bai_older_than_data.bam"

# When using the index, samtools outputs a warning
validate $(samtools idxstats "$DIR_OUT/bai_older_than_data.bam" 2>&1 | grep "The index file is older than the data file")
