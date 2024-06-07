#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good FASTA files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/fasta/good
DIR_BASIC="$DIR_OUT/basic_dna.fa"

# ------------------------------------------------------------------------------
# Compressed. Use bgzip because gzip adds timestamp in file, which makes git show
# diffs each time this script runs.
# ------------------------------------------------------------------------------

log "Creating compressed FASTA file"
cp "$DIR_BASIC" "$DIR_OUT/compressed.fa"
bgzip -f "$DIR_OUT/compressed.fa"
validate "$(diff "$DIR_BASIC" <(gunzip -c "$DIR_OUT/compressed.fa.gz") && echo "ok")"

# ------------------------------------------------------------------------------
# Multiline FASTA
# ------------------------------------------------------------------------------

log "Creating FASTA split across multiple lines"
seqtk seq -l 20 "$DIR_BASIC" > "$DIR_OUT/multiline.fa"
validate "$(diff <(wc -l < "$DIR_OUT/multiline.fa") <(wc -l < "$DIR_BASIC"))"

# ------------------------------------------------------------------------------
# Sequence names with spaces
# ------------------------------------------------------------------------------

log "Creating FASTA with spaces in sequence name"
seqtk rename "$DIR_BASIC" "prefix  with  spaces" > "$DIR_OUT/name_contains_spaces.fa"
validate "$(diff <(grep -c " " "$DIR_BASIC") <(grep -c " " "$DIR_OUT/name_contains_spaces.fa"))"

# ------------------------------------------------------------------------------
# Duplicate sequence names
# ------------------------------------------------------------------------------

log "Creating FASTA with duplicate sequence names"
sed 's/sequence1/sequence2/g' "$DIR_BASIC" > "$DIR_OUT/duplicate_sequence_names.fa"
validate "$(diff <(grep '>' "$DIR_BASIC") <(grep '>' "$DIR_OUT/duplicate_sequence_names.fa"))"

# ------------------------------------------------------------------------------
# Spaces between sequences...
# ------------------------------------------------------------------------------

log "Creating FASTA with empty lines between sequences"
sed 's/>sequence/\n>sequence/g' "$DIR_BASIC" > "$DIR_OUT/empty_lines.fa"
validate "$(diff <(wc -l < "$DIR_BASIC") <(wc -l < "$DIR_OUT/empty_lines.fa"))"
