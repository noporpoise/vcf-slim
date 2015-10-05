
# You need to pass the path to htslib e.g.:
#   make HTSLIB=~/bioinf/htslib/

HTSLIB=../htslib

CFLAGS=-Wall -Wextra -O2
LIBS=-lz -lm -lpthread

all: bin/vcfhp bin/vcfdist bin/vcfcontigs

bin/%: src/%.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ -I $(HTSLIB)/htslib $< $(HTSLIB)/libhts.a $(LIBS)

$(DIRS):
	mkdir -p bin

clean:
	rm -rf bin

test: bin/vcfhp
	rm -rf tests/ref.fa.fai
	bin/vcfhp tests/ref.fa tests/in.vcf > tests/out.vcf
	diff -q tests/out.vcf tests/ans.vcf
	@echo Looks good.

.PHONY: all clean test
