
# You need to pass the path to htslib e.g.:
#   make HTSLIB=~/bioinf/htslib/

HTSLIB=../htslib

CFLAGS=-Wall -Wextra -O2
LIBS=-lz -lm

all: vcfhp

%: %.c
	$(CC) $(CFLAGS) -I $(HTSLIB)/htslib -o $@ $< $(HTSLIB)/libhts.a $(LIBS)

clean:
	rm -rf vcfhp

test: vcfhp
	rm -rf tests/ref.fa.fai
	./vcfhp tests/ref.fa tests/in.vcf > tests/out.vcf
	diff -q tests/out.vcf tests/ans.vcf
	@echo Looks good.

.PHONY: all clean test
