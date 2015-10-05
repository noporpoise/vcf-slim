#vcfhp

Annotate a VCF with homopolymer run tag (`HRun=`)

`HRun` is the number of bases either side of the variant in the reference that
match the variant sequence. Annotation is only applied to biallelic indels.

Compile:

    make HTSLIB=../htslib all

Run tests:

    make test

Run:

    bin/vcfhp ref.fa in.vcf > out.vcf

Running will generate an index for `ref.fa` (at `ref.fa.fai`) if it does not
already exist. Adds info tag e.g. `HRun=4` to biallelic indels on homopolyer
runs greater than one.

Useful page on Homopolymer Run annotation from joinx:
https://github.com/genome/joinx/tree/42a7caf2c75923ccbbfdcd072673586cda684cac/integration-test/data/vcf-annotate-homopolymers
