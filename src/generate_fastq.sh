#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good FASTQ files
# ==============================================================================

URL_R1="https://42basepairs.com/download/r2/genomics-data/reads_NA12878_R1.fastq.gz"
URL_R2="https://42basepairs.com/download/r2/genomics-data/reads_NA12878_R2.fastq.gz"
URL2_R1="https://42basepairs.com/download/s3/1000genomes/technical/working/20150612_chm1_data/sequence_read/SRR642626_1.filt.fastq.gz"
URL2_R2="https://42basepairs.com/download/s3/1000genomes/technical/working/20150612_chm1_data/sequence_read/SRR642626_2.filt.fastq.gz"
DIR_OUT=$DIR_SRC/../data/fastq/good
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Basic
# ------------------------------------------------------------------------------

log "Creating valid FASTQ file"
DIR_BASIC="$DIR_OUT/basic.fastq"
curl -s "$URL_R1" | seqtk seq -l 0 | head -n 12 > "$DIR_BASIC"
validate "ok"

# ------------------------------------------------------------------------------
# Quality line starts with "@", which can trip up FASTQ parsers
# ------------------------------------------------------------------------------

log "Creating FASTQ file where the quality line (line 4) starts with @"
sed '4s/./@/' "$DIR_BASIC" > "$DIR_OUT/quality_@.fastq"
validate "$(diff <(echo 3) <(grep -c "^@" "$DIR_OUT/quality_@.fastq"))"

# ------------------------------------------------------------------------------
# Sequence/quality on multiple lines, which breaks parsers that assume 4 lines per read
# ------------------------------------------------------------------------------

log "Creating FASTQ split across multiple lines"
seqtk seq -l 60 "$DIR_BASIC" > "$DIR_OUT/multiline.fastq"
validate "$(diff <(wc -l < "$DIR_OUT/multiline.fastq") <(wc -l < "$DIR_BASIC"))"

# ------------------------------------------------------------------------------
# Non-empty "+" line
# ------------------------------------------------------------------------------

log "Creating FASTQ with a + line with duplicate read names"
awk 'BEGIN { read = -1 } { if((NR-1) % 4 == 0) { read=$0; } if($0 == "+") { sub("@", "+", read); print(read) } else { print }  }' "$DIR_BASIC" > "$DIR_OUT/duplicate_+.fastq"
validate "$(diff "$DIR_BASIC" "$DIR_OUT/duplicate_+.fastq" )"

# ------------------------------------------------------------------------------
# Interleaved R1 and R2 reads
# ------------------------------------------------------------------------------

log "Creating interleaved FASTQ (read name does not contain R1/R2)"
seqtk mergepe "$DIR_BASIC" <(curl -s "$URL_R2" | seqtk seq -l 0 | head -n 12) > "$DIR_OUT/interleaved.fastq"
validate "$(diff <(wc -l < "$DIR_BASIC") <(wc -l "$DIR_OUT/interleaved.fastq"))"

log "Creating interleaved FASTQ (read name contains R1/R2)"
seqtk mergepe <(curl -s "$URL2_R1" | seqtk seq -l 0 | head -n 12) <(curl -s "$URL2_R2" | seqtk seq -l 0 | head -n 12) > "$DIR_OUT/interleaved2.fastq"
validate "$(diff <(wc -l < "$DIR_BASIC") <(wc -l "$DIR_OUT/interleaved2.fastq"))"


# ==============================================================================
# Bad FASTQ files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/fastq/bad
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Truncated
# ------------------------------------------------------------------------------

log "Creating truncated FASTQ (missing lines)"
head -n 10 "$DIR_BASIC" > "$DIR_OUT/truncated_clean.fastq"
validate "$(diff <(seqtk seq "$DIR_BASIC" | wc -l) <(seqtk seq "$DIR_OUT/truncated_clean.fastq" | wc -l))"

log "Creating truncated FASTQ (half-way through a line)"
head -c 200 "$DIR_BASIC" > "$DIR_OUT/truncated_halfway.fastq"
validate "$(diff <(seqtk seq "$DIR_BASIC" | wc -l) <(seqtk seq "$DIR_OUT/truncated_halfway.fastq" | wc -l))"

# ------------------------------------------------------------------------------
# Not truncated but a read in the middle has len(sequence) != len(quality)
# ------------------------------------------------------------------------------

log "Creating FASTQ file where we're missing quality scores for some bases"
sed '8s/.....//' "$DIR_BASIC" > "$DIR_OUT/quality_mismatch.fastq"
validate "$(diff <(seqtk seq "$DIR_BASIC" | wc -l) <(seqtk seq "$DIR_OUT/truncated_clean.fastq" | wc -l))"
