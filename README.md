Annotate a VCF with homopolymer run tag (`HRun=`)

Compile:

    make HTSLIB=../htslib

Run:

    ./vcfhp ref.fa in.vcf > out.vcf

Running will generate an index for `ref.fa` (at `ref.fa.fai`) if it does not
already exist.
