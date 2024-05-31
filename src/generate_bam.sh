#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good BAM files
# ==============================================================================

URL="https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
DIR_OUT=$DIR_SRC/../data/bam/good
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Basic BAM file
# ------------------------------------------------------------------------------

log "Creating valid BAM file"
samtools view -b -h "$URL" 11:82365011-82366010 > "$DIR_OUT/basic.bam"
validate "$(samtools quickcheck -vvv "$DIR_OUT/basic.bam" 2>&1 | grep "good EOF block")"

# ------------------------------------------------------------------------------
# Unsorted BAM file
# ------------------------------------------------------------------------------

log "Creating unsorted BAM file"
samtools view -b -h "$URL" 11:128989445-128990444 11:82365011-82366010 > "$DIR_OUT/unsorted.bam"
validate "$(samtools index "$DIR_OUT/unsorted.bam" 2>&1 | grep "Unsorted positions on sequence")"

# ------------------------------------------------------------------------------
# No mapped data in BAM
# ------------------------------------------------------------------------------

log "Creating BAM with no mapped reads"
samtools reset "$DIR_OUT/basic.bam" -O BAM > "$DIR_OUT/no_mapped_reads.bam"
validate "$(samtools quickcheck "$DIR_OUT/no_mapped_reads.bam" 2>&1 | grep "no targets in header")"


# ==============================================================================
# Bad BAM files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/bam/bad
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Truncated BAM file
# ------------------------------------------------------------------------------

log "Creating truncated BAM file"
curl -s -H "Range: bytes=0-25000" "$URL" > "$DIR_OUT/truncated.bam"
validate "$(samtools view "$DIR_OUT/truncated.bam" 2>&1 | grep "EOF marker is absent. The input is probably truncated")"

# ------------------------------------------------------------------------------
# BAI index is older than the BAM file (might needs reindexing to reflect latest BAM file contents)
# ------------------------------------------------------------------------------

log "Creating BAI index that is older than the BAM file"
samtools view -b -h "$URL" 11:82365011-82366010 > "$DIR_OUT/bai_older_than_data.bam"
samtools index "$DIR_OUT/bai_older_than_data.bam"
sleep 1  # wait before run touch to make sure the timestamp is different enough
touch "$DIR_OUT/bai_older_than_data.bam"
validate "$(samtools idxstats "$DIR_OUT/bai_older_than_data.bam" 2>&1 | grep "The index file is older than the data file")"
