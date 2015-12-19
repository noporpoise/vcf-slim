#!/bin/bash

set -eou pipefail

DATE=$(date '+%Y%m%d')

echo '##fileformat=VCFv4.1
##fileDate='"$DATE"'
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'

echo '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT' | tr ' ' "\t"

# Input from vcf-slim/bin/vcfcontigs:
# NC_009648.1:813:.:C:T:1:pos:50:ltrim:0:rtrim:0
gzip -fcd $@ | tr ':' '\t' | cut -f 1-5 | \
  while read line; do printf "$line\t.\t.\tGT\n"; done
