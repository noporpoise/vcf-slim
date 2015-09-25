#vcfhp

Annotate a VCF with homopolymer run tag (`HRun=`)

Homopolymer is defined as a variant which only changes the length of a run of bases in the reference. 

Compile:

    make HTSLIB=../htslib

Run tests:

    make test

Run:

    ./vcfhp ref.fa in.vcf > out.vcf

Running will generate an index for `ref.fa` (at `ref.fa.fai`) if it does not
already exist. Adds info tag e.g. `HRun=4` to biallelic indels on homopolyer
runs greater than one.

Useful page on Homopolymer Run annotation from joinx:
https://github.com/genome/joinx/tree/42a7caf2c75923ccbbfdcd072673586cda684cac/integration-test/data/vcf-annotate-homopolymers

##TODO:

Currently only deals with single base repeats. Multi-base repeats should be
added.
