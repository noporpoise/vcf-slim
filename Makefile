
# You need to pass the path to htslib e.g.:
#   make HTSLIB=~/bioinf/htslib/

HTSLIB=../htslib
LIBS=-lz -lm -lpthread
CFLAGS=-Wall -Wextra $(OPTS)

ifdef DEBUG
	OPTS=-O0 -g
else
	OPTS=-O2
endif

all: bin/vcfhp bin/vcfdist bin/vcfcontigs

bin/%: src/%.c src/common.h
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
