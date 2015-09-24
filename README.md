Annotate a VCF with homopolymer run tag (`HRun=`)

Homopolymer is defined as a variant which only changes the length of a run of bases in the reference. 

Compile:

    make HTSLIB=../htslib

Run tests:

    make test

Run:

    ./vcfhp ref.fa in.vcf > out.vcf

Running will generate an index for `ref.fa` (at `ref.fa.fai`) if it does not
already exist.
