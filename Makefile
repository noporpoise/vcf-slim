
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

.PHONY: all clean
