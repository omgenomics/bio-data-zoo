#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good BED files
# ==============================================================================

URL="https://42basepairs.com/download/r2/genomics-data/regions_CHM13.bed.gz"
DIR_OUT=$DIR_SRC/../data/bed/good
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Basic
# ------------------------------------------------------------------------------

log "Creating valid BED file"
DIR_BASIC="$DIR_OUT/basic.bed"
curl -s "$URL" | zcat | head -n 10 > "$DIR_BASIC"
validate "ok"

# ------------------------------------------------------------------------------
# Unsorted intervals
# ------------------------------------------------------------------------------

log "Creating unsorted BED file"
curl -s "$URL" | zcat | head -n 50 | tail -n 10 > "$DIR_OUT/unsorted.bed"
cat "$DIR_BASIC" >> "$DIR_OUT/unsorted.bed"
validate "$(diff "$DIR_OUT/unsorted.bed" <(bedtools sort -i "$DIR_OUT/unsorted.bed"))"

# ------------------------------------------------------------------------------
# Compressed
# ------------------------------------------------------------------------------
# Note that tabix requires bgzip compression. Also, using `gzip` adds a timestamp
# to the header, which makes it look different to git every time the script runs.
# ------------------------------------------------------------------------------

log "Creating compressed BED file (bgzip)"
cp "$DIR_BASIC" "$DIR_OUT/compressed.bed"
bgzip -f "$DIR_OUT/compressed.bed"
validate "ok"

log "Creating indexed BED.bgz file (TBI, CSI)"
cp "$DIR_BASIC" "$DIR_OUT/indexed_tbi.bed"
cp "$DIR_BASIC" "$DIR_OUT/indexed_csi.bed"
bgzip -f "$DIR_OUT/indexed_tbi.bed"
bgzip -f "$DIR_OUT/indexed_csi.bed"
tabix -p bed "$DIR_OUT/indexed_tbi.bed.gz"
tabix -p bed --csi "$DIR_OUT/indexed_csi.bed.gz"
validate "ok"


# ==============================================================================
# Bad BED files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/bed/bad
mkdir -p "$DIR_OUT"

# ------------------------------------------------------------------------------
# Spaces instead of tabs
# ------------------------------------------------------------------------------

log "Creating BED file using spaces instead of tabs"
sed 's/\t/    /g' "$DIR_BASIC" > "$DIR_OUT/spaces.bed"
validate "$(bedtools merge -i "$DIR_OUT/spaces.bed" 2>&1 | grep "unable to open file or unable to determine types")"

# ------------------------------------------------------------------------------
# Invalid ranges
# ------------------------------------------------------------------------------

log "Creating BED file with negative coordinates"
sed 's/3634/-3634/' "$DIR_BASIC" > "$DIR_OUT/negative_coords.bed"
validate "$(bedtools merge -i "$DIR_OUT/negative_coords.bed" 2>&1 | grep "Invalid record in file")"

log "Creating BED file with start > end coordinates"
sed 's/3634/9999/' "$DIR_BASIC" > "$DIR_OUT/start_greater_than_end_coords.bed"
validate "$(bedtools merge -i "$DIR_OUT/start_greater_than_end_coords.bed" 2>&1 | grep "unable to open file or unable to determine types")"

log "Creating BED file with non-integer coordinates"
sed 's/3634/3.63/' "$DIR_BASIC" > "$DIR_OUT/non_integer_coords.bed"
validate "$(bedtools merge -i "$DIR_OUT/non_integer_coords.bed" 2>&1 | grep "unable to open file or unable to determine types")"
