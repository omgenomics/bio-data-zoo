#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good FASTQ files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/fastq/good
DIR_BASIC_R1="$DIR_OUT/basic_R1.fastq"
DIR_BASIC_R2="$DIR_OUT/basic_R2.fastq"

# ------------------------------------------------------------------------------
# Compressed
# ------------------------------------------------------------------------------

log "Creating compressed FASTQ file"
cp "$DIR_BASIC_R1" "$DIR_OUT/compressed.fastq"
bgzip -f "$DIR_OUT/compressed.fastq"
validate "$(diff "$DIR_BASIC_R1" <(gunzip -c "$DIR_OUT/compressed.fastq") && echo "ok")"

# ------------------------------------------------------------------------------
# Quality line starts with "@", which can trip up FASTQ parsers
# ------------------------------------------------------------------------------

log "Creating FASTQ file where the quality line (line 4) starts with @"
sed '4s/./@/' "$DIR_BASIC_R1" > "$DIR_OUT/quality_@.fastq"
validate "$(diff <(echo 3) <(grep -c "^@" "$DIR_OUT/quality_@.fastq"))"

# ------------------------------------------------------------------------------
# Sequence/quality on multiple lines, which breaks parsers that assume 4 lines per read.
# This is quite rare but is technically a valid FASTQ file.
# ------------------------------------------------------------------------------

log "Creating FASTQ split across multiple lines"
seqtk seq -l 20 "$DIR_BASIC_R1" > "$DIR_OUT/multiline.fastq"
validate "$(diff <(wc -l < "$DIR_OUT/multiline.fastq") <(wc -l < "$DIR_BASIC_R1"))"

# ------------------------------------------------------------------------------
# Non-empty "+" line
# ------------------------------------------------------------------------------

log "Creating FASTQ with a + line with duplicate read names"
awk 'BEGIN { read = -1 } { if((NR-1) % 4 == 0) { read=$0; } if($0 == "+") { sub("@", "+", read); print(read) } else { print }  }' "$DIR_BASIC_R1" > "$DIR_OUT/duplicate_+.fastq"
validate "$(diff "$DIR_BASIC_R1" "$DIR_OUT/duplicate_+.fastq" )"

# ------------------------------------------------------------------------------
# Interleaved R1 and R2 reads
# ------------------------------------------------------------------------------

log "Creating interleaved FASTQ (read name does not contain R1/R2)"
seqtk mergepe "$DIR_BASIC_R1" "$DIR_BASIC_R2" > "$DIR_OUT/interleaved.fastq"
validate "$(diff <(wc -l < "$DIR_BASIC_R1") <(wc -l "$DIR_OUT/interleaved.fastq"))"


# ==============================================================================
# Bad FASTQ files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/fastq/bad
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Truncated
# ------------------------------------------------------------------------------

log "Creating truncated FASTQ (missing lines)"
head -n 10 "$DIR_BASIC_R1" > "$DIR_OUT/truncated_clean.fastq"
validate "$(diff <(seqtk seq "$DIR_BASIC_R1" | wc -l) <(seqtk seq "$DIR_OUT/truncated_clean.fastq" | wc -l))"

log "Creating truncated FASTQ (half-way through a line)"
head -c 200 "$DIR_BASIC_R1" > "$DIR_OUT/truncated_halfway.fastq"
validate "$(diff <(seqtk seq "$DIR_BASIC_R1" | wc -l) <(seqtk seq "$DIR_OUT/truncated_halfway.fastq" | wc -l))"

# ------------------------------------------------------------------------------
# Not truncated but a read in the middle has len(sequence) != len(quality)
# ------------------------------------------------------------------------------

log "Creating FASTQ file where we're missing quality scores for some bases"
sed '8s/.....//' "$DIR_BASIC_R1" > "$DIR_OUT/quality_mismatch.fastq"
validate "$(diff <(seqtk seq "$DIR_BASIC_R1" | wc -l) <(seqtk seq "$DIR_OUT/truncated_clean.fastq" | wc -l))"
