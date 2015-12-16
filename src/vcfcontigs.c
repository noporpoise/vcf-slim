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
"\n"
"  Contigs are named:\n"
"    CHROM:POS:ID:REF:ALT:ALTIDX:pos:OFFSET:ltrim:LTRIM:rtrim:RTRIM\n"
"\n"
"  Options:\n"
"   -t,--trim         Remove matching bases from left and right of alleles\n"
"   -R,--no-ref       Don't print ref contig\n"
"   -A,--max-alt <A>  Skip ALTs longer than <A>\n";

static void print_usage()
{
  fprintf(stderr, "usage: %s [options] <flank> <ref.fa> <in.vcf>\n", cmd);
  fprintf(stderr, usage);
  exit(-1);
}

static const char shortopts[] = "tRA:";

static struct option longopts[] =
{
  {"trim",          no_argument, NULL, 't'},
  {"no-ref",        no_argument, NULL, 'R'},
  {"max-alt", required_argument, NULL, 'A'},
  {NULL, 0, NULL, 0}
};

size_t n_max_ref_skipped = 0, n_max_alt_skipped = 0;

static inline void print_var(bcf_hdr_t *hdr, bcf1_t *v,
                             int flank, bool trim_seqs, int max_alt,
                             bool print_ref, bool print_alt,
                             char *chrom, int chromlen)
{
  assert(chrom);

  FILE *fout = stdout;
  char *str, *estr, *astr;
  int start, end, j;
  size_t ltrim = 0, rtrim = 0, trim, pos, rlen, alen, i;
  size_t first = print_ref ? 0 : 1, last = print_alt ? v->n_allele : 1;

  // trim alleles
  if(trim_seqs) trim_alleles(v, &ltrim, &rtrim);
  trim = ltrim + rtrim;
  pos = v->pos + ltrim;
  // fprintf(stderr, "trim:%zu rlen:%i\n", trim, v->rlen);
  assert((size_t)v->rlen >= trim);
  rlen = v->rlen - trim;
  start = MAX2(0, pos - flank);
  end = MIN2(chromlen, (int)(pos + rlen + flank));

  if(max_alt > 0) {
    if(rlen > (size_t)max_alt) {
      n_max_ref_skipped++;
      return;
    }
    for(j = 1; j < v->n_allele; j++)
      if(strlen(v->d.allele[j])-trim <= (size_t)max_alt)
        break;
    if(j == v->n_allele) {
      n_max_alt_skipped += v->n_allele-1;
      return;
    }
  }

  for(i = first; i < last; i++)
  {
    astr = v->d.allele[i];
    alen = strlen(v->d.allele[i]);
    if(max_alt > 0 && alen > (size_t)max_alt) {
      n_max_alt_skipped++;
      continue;
    }
    fprintf(fout, ">%s:%i:%s:%s:%s:%zu:pos:%zu:ltrim:%zu:rtrim:%zu\n",
            bcf_seqname(hdr, v), v->pos+1, v->d.id, v->d.allele[0], v->d.allele[i],
            i, pos-start, ltrim, rtrim);
    for(str = chrom+start, estr = chrom+pos; str < estr; str++) fputc(*str, fout);
    for(str = astr+ltrim, estr = astr+alen-trim; str < estr; str++) fputc(*str, fout);
    for(str = chrom+pos+rlen, estr = chrom+end; str < estr; str++) fputc(*str, fout);
    fputc('\n', fout);
  }
}

int main(int argc, char **argv)
{
  cmd = argv[0];
  bool trim_alleles = false, print_ref = true, print_alt = true;
  int max_alt = -1;

  // Arg parsing
  int c;
  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(c) {
      case 0: /* flag set */ break;
      case 't': trim_alleles = true; break;
      case 'R': print_ref = false; break;
      case 'A': max_alt = atoi(optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        die("Bad option: %s", argv[optind-1]);
      default:
        die("Coder error");
    }
  }

  if(optind + 3 != argc) print_usage();

  char **args = argv + optind;
  int flank = atoi(args[0]);
  const char *ref_path = args[1];
  const char *vcf_path = args[2];
  const char *out_path = "-";

  fprintf(stderr, "[vcfcontigs] flank: %i\n", flank);
  fprintf(stderr, "[vcfcontigs] ref: %s\n", ref_path);
  fprintf(stderr, "[vcfcontigs] vcf: %s\n", vcf_path);
  fprintf(stderr, "[vcfcontigs] out: %s\n", out_path);
  fprintf(stderr, "[vcfcontigs] max-alt: %i\n", max_alt);

  // open input / output
  htsFile *vin = hts_open(vcf_path, "r");
  if(vin == NULL) die("Cannot read %s: %s", vcf_path, strerror(errno));
  // htsFile *vout = hts_open(out_path, "w");
  // if(vout == NULL) die("Cannot write %s: %s", out_path, strerror(errno));

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
    print_var(hdr, v, flank, trim_alleles, max_alt,
              print_ref, print_alt,
              chrom, chromlen);
  }

  free(chrom);
  fai_destroy(fai);

  bcf_hdr_destroy(hdr);
  hts_close(vin);
  // hts_close(vout);

  fprintf(stderr, "[vcfcontigs] skipped %zu refs\n", n_max_ref_skipped);
  fprintf(stderr, "[vcfcontigs] skipped %zu alts\n", n_max_alt_skipped);

  return EXIT_SUCCESS;
}
