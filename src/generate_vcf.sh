#!/bin/bash

DIR_SRC="${BASH_SOURCE%/*}"
source "$DIR_SRC/lib.sh" || exit

# ==============================================================================
# Good VCF files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/vcf/good
DIR_BASIC="$DIR_OUT/basic.vcf"

# Compressed
log "Creating compressed VCF file (bgzip)"
cp "$DIR_BASIC" "$DIR_OUT/compressed.vcf"
bgzip -f "$DIR_OUT/compressed.vcf"
validate "$(diff "$DIR_BASIC" <(gunzip -c "$DIR_OUT/compressed.vcf.gz") && echo "ok")"

# Indexed VCF
log "Creating indexed VCF file (TBI, CSI)"
cp "$DIR_OUT/compressed.vcf.gz" "$DIR_OUT/indexed_csi.vcf.gz"
cp "$DIR_OUT/compressed.vcf.gz" "$DIR_OUT/indexed_tbi.vcf.gz"
bcftools index --csi "$DIR_OUT/indexed_csi.vcf.gz"
bcftools index --tbi "$DIR_OUT/indexed_tbi.vcf.gz"
validate "ok"

# Indexed BCF
log "Creating indexed BCF file (CSI)"
cp "$DIR_OUT/basic.bcf" "$DIR_OUT/indexed.bcf"
bcftools index "$DIR_OUT/indexed.bcf"
validate "ok"

# ==============================================================================
# Bad VCF files
# ==============================================================================

DIR_OUT=$DIR_SRC/../data/vcf/bad

# Missing an INFO field (bcftools outputs warnings)
log "Creating VCF with a missing INFO field"
sed "/^##INFO=<ID=AN,/d;" "$DIR_BASIC" > "$DIR_OUT/missing_info_field.vcf"
validate "$(bcftools view "$DIR_OUT/missing_info_field.vcf" 2>&1 | grep "W::vcf_parse_info")"
