#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good FASTA files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/fasta/good
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Basic
# ------------------------------------------------------------------------------

log "Creating valid FASTA file"
DIR_BASIC="$DIR_OUT/basic.fa"
cat > "$DIR_BASIC" <<EOF
>sequence1
TTAACTCTTGCTGTCCCNNNNNCAGCAGAGGCGCCAATCATAGTTTTAAGAGAAGGAGGAAGGGACTGAGTGAGGGGAGAGTGGATGAGTAATTAGGTGG
>sequence2
CAAGTGCAGTTGCCTGGATTTTGCATTGCAAGATTCCCAGGGTGTACTTTTCACATGCATTACAACGTGAATGATCTGCCCTACCATGGTTGTGTAGTAT
>sequence3
CCAGGGGCACCCATGCTGGACCAGGCTGCACGGCCAGGGCTGCTGGGACAGAGGTAGGTGGGGTCCTGGAGCAGAGCCAGCTCCACACGCCAATAGCCCA
EOF
validate "ok"

# ------------------------------------------------------------------------------
# Basic (protein)
# ------------------------------------------------------------------------------

log "Creating valid FASTA file (protein)"
cat > "$DIR_OUT/protein.fa" <<EOF
>sequence1
LYLIFGAWAGMVGTALSLLIRAELGQPGTLLGDDQIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMI
>sequence2
GAPDMAFPRMNNMSFWLLPPSFLLLLASSTVEAGAGTGWTVYPPLAGNLAHAGASVDLAIFSLHLAGVSS
>sequence3
ILGAINFITTAINMKPPTLSQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGD
EOF
validate "ok"

# ------------------------------------------------------------------------------
# Basic (aligned)
# ------------------------------------------------------------------------------

log "Creating aligned FASTA file (- for gaps)"
cat > "$DIR_OUT/aligned.fa" <<EOF
>sequence1
ATCTACGATCGAGCTACT
>sequence2
ATC----ATCGACCCACT
EOF
validate "ok"

# ------------------------------------------------------------------------------
# Compressed. Use bgzip because gzip adds timestamp in file, which makes git show
# diffs each time this script runs.
# ------------------------------------------------------------------------------

log "Creating compressed FASTA file"
cp "$DIR_BASIC" "$DIR_OUT/compressed.fa"
bgzip -f "$DIR_OUT/compressed.fa"
validate "ok"

# ------------------------------------------------------------------------------
# Multiline FASTA
# ------------------------------------------------------------------------------

log "Creating FASTA split across multiple lines"
seqtk seq -l 60 "$DIR_BASIC" > "$DIR_OUT/multiline.fa"
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
