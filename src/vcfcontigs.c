#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <assert.h>
#include <getopt.h>

#include "common.h"

#include "hts.h"
#include "vcf.h"
#include "faidx.h"

const char *cmd = NULL;

const char usage[] = 
"  Print <flank> bp either side of a variant, as FASTA to STDOUT.\n"
"  Options:\n"
"   -t,--trim      Remove matching bases from left and right of alleles\n"
"   -R,--no-ref    Don't print ref contig\n";

static void print_usage()
{
  fprintf(stderr, "usage: %s [options] <flank> <ref.fa> <in.vcf>\n", cmd);
  fprintf(stderr, usage);
  exit(-1);
}

static const char shortopts[] = "tR";

static struct option longopts[] =
{
  {"trim",    no_argument, NULL, 't'},
  {"no-ref",  no_argument, NULL, 'R'},
  {NULL, 0, NULL, 0}
};

static inline void print_var(bcf_hdr_t *hdr, bcf1_t *v,
                             int flank, bool print_ref, bool print_alt,
                             bool trim_seqs,
                             char *chrom, int chromlen)
{
  assert(chrom);

  FILE *fout = stdout;
  char *str, *estr, *astr;
  int start, end;
  size_t ltrim = 0, rtrim = 0, trim, pos, rlen, alen, i;
  size_t first = print_ref ? 0 : 1, last = print_alt ? v->n_allele : 1;

  // trim alleles
  if(trim_seqs) trim_alleles(v, &ltrim, &rtrim);
  trim = ltrim + rtrim;
  pos = v->pos + ltrim;
  rlen = v->rlen - trim;
  start = MAX2(0, pos - flank);
  end = MIN2(chromlen, (int)(pos + rlen + flank));

  for(i = first; i < last; i++)
  {
    astr = v->d.allele[i];
    alen = strlen(v->d.allele[i]);
    fprintf(fout, ">%s_%i_%s_%s_%s_%zu_%zu_%zu\n", bcf_seqname(hdr, v), v->pos,
            v->d.id, v->d.allele[0], v->d.allele[i], i, rtrim, ltrim);
    for(str = chrom+start, estr = chrom+pos; str < estr; str++) fputc(*str, fout);
    for(str = astr+ltrim, estr = astr+alen-trim; str < estr; str++) fputc(*str, fout);
    for(str = chrom+pos+rlen, estr = chrom+end; str < estr; str++) fputc(*str, fout);
    fputc('\n', fout);
  }
}

int main(int argc, char **argv)
{
  bool trim_alleles = false, print_ref = true, print_alt = true;

  // Arg parsing
  int c;
  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(c) {
      case 0: /* flag set */ break;
      case 't': trim_alleles = true; break;
      case 'R': print_ref = false; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("Bad option: %s", argv[optind-1]);
      default:
        die("Coder error");
    }
  }

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
    print_var(hdr, v, flank, print_ref, print_alt, trim_alleles,
              chrom, chromlen);
  }

  free(chrom);
  fai_destroy(fai);

  bcf_hdr_destroy(hdr);
  hts_close(vin);
  hts_close(vout);

  return EXIT_SUCCESS;
}
