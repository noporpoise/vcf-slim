#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <assert.h>

#include "hts.h"
#include "vcf.h"
#include "faidx.h"

const char *cmd = NULL;

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0)

static void print_usage()
{
  fprintf(stderr, "usage: %s <flank> <ref.fa> <in.vcf>\n", cmd);
  fprintf(stderr, "  Print <flank> bp either side of a variant\n");
  fprintf(stderr, "  Print FASTA to STDOUT.\n");
  exit(-1);
}

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MAX2(a,b) ((a) >= (b) ? (a) : (b))

static inline void fetch_chrom(const bcf_hdr_t *hdr, bcf1_t *v,
                               faidx_t *fai, int *refid,
                               char **chrom, int *chromlen)
{
  if(*refid != v->rid) {
    free(*chrom);
    *chrom = fai_fetch(fai, bcf_seqname(hdr, v), chromlen);
    if(*chrom == NULL) die("Cannot find chrom '%s'", bcf_seqname(hdr, v));
    *refid = v->rid;
  }
}

static inline void print_var(bcf_hdr_t *hdr, bcf1_t *v, int flank,
                             char *chrom, int chromlen)
{
  assert(chrom);

  FILE *fout = stdout;
  char *str, *estr;
  int i, start, end;

  start = MAX2(0, v->pos - flank);
  end = MIN2(chromlen, v->pos + v->rlen + flank);

  for(i = 0; i < v->n_allele; i++) {
    fprintf(fout, ">%s_%s_%i_%i\n", v->d.id, bcf_seqname(hdr, v), v->pos, i);
    for(str = chrom+start, estr = chrom+v->pos; str < estr; str++)
      fputc(*str, fout);
    fputs(v->d.allele[i], fout);
    for(str = chrom+v->pos+v->rlen, estr = chrom+end; str < estr; str++)
      fputc(*str, fout);
    fputc('\n', fout);
  }
}

int main(int argc, char **argv)
{
  cmd = argv[0];
  if(argc != 4) print_usage();

  int flank = atoi(argv[1]);
  const char *ref_path = argv[2];
  const char *vcf_path = argv[3];
  const char *out_path = "-";

  // open input / output
  htsFile *vin = hts_open(vcf_path, "r");
  if(vin == NULL) die("Cannot read %s: %s", vcf_path, strerror(errno));
  htsFile *vout = hts_open(out_path, "w");
  if(vout == NULL) die("Cannot write %s: %s", out_path, strerror(errno));

  // read header
  bcf_hdr_t *hdr = bcf_hdr_read(vin);

  // open ref
  faidx_t *fai = fai_load(ref_path);
  if(fai == NULL) die("Build %s.fai: samtools faidx %s", ref_path, ref_path);

  // Load ref chrom
  int refid = -1, chromlen = 0;
  char *chrom = NULL;

  // read vcf entries
  bcf1_t *v = bcf_init();

  while(bcf_read(vin, hdr, v) >= 0)
  {
    bcf_unpack(v, BCF_UN_STR);
    fetch_chrom(hdr, v, fai, &refid, &chrom, &chromlen);
    print_var(hdr, v, flank, chrom, chromlen);
  }

  free(chrom);
  fai_destroy(fai);

  bcf_hdr_destroy(hdr);
  hts_close(vin);
  hts_close(vout);

  return EXIT_SUCCESS;
}
