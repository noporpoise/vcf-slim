#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

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
  fprintf(stderr, "usage: %s <ref.fa> <in.vcf>\n", cmd);
  fprintf(stderr, "  Add homopolymer run annotations to VCF files (HRun=).\n");
  fprintf(stderr, "  Prints to STDOUT.\n");
  exit(-1);
}

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MAX2(a,b) ((a) >= (b) ? (a) : (b))

// Call bases_match() and check reflen != altlen
// Note: hp may be less than variant length
// pos is relative start in reg
static inline int calc_homopolymer(char *ref, int reflen, char *alt, int altlen,
                                   int pos, char *reg, int reglen)
{
  int i, k, hp = 0;
  const char *r = reflen ? ref : alt;
  int rl = reflen ? reflen : altlen;
  for(i = pos-1, k=rl-1; i >= 0 && reg[i] == r[k]; i--, k = k ? k-1 : rl-1) hp++;
  for(i = pos+reflen, k=0; i < reglen && reg[i] == r[k]; i++, k++, k=k<rl?k:0) hp++;
  return hp;
}

int main(int argc, char **argv)
{
  cmd = argv[0];
  if(argc != 3) print_usage();
  const char *ref_path = argv[1];
  const char *vcf_path = argv[2];

  const char *out_path = "-";

  // open input / output
  htsFile *vin = hts_open(vcf_path, "r");
  if(vin == NULL) die("Cannot read %s: %s", vcf_path, strerror(errno));
  htsFile *vout = hts_open(out_path, "w");
  if(vout == NULL) die("Cannot write %s: %s", out_path, strerror(errno));

  // Note: Build ref.fa.fai manually with: samtools faidx ref.fa

  // open ref
  faidx_t *fai = fai_load(ref_path);
  if(fai == NULL) die("Build %s.fai: samtools faidx %s", ref_path, ref_path);

  // read/write headers
  bcf_hdr_t *hdrin = bcf_hdr_read(vin);
  bcf_hdr_t *hdrout = bcf_hdr_dup(hdrin);
  bcf_hdr_append(hdrout, "##INFO=<ID=HRun,Number=1,Type=Integer,Description=\""
                         "Homopolymer run in ref in bp (not including variant)\">\n");
  if(bcf_hdr_write(vout, hdrout) != 0) die("Cannot write header");

  // read vcf entries
  bcf1_t *v = bcf_init();
  int hp, window = 100;

  while(bcf_read(vin, hdrin, v) >= 0)
  {
    // Unpack ref,alt,filter,info info
    bcf_unpack(v, BCF_UN_SHR);

    // Only annotate if biallelic
    if(v->n_allele == 2)
    {
      char *ref = v->d.allele[0];
      char *alt = v->d.allele[1];

      int reflen = strlen(ref);
      int altlen = strlen(alt);
      int pos = v->pos;

      // Left trim
      while(reflen && altlen && *ref == *alt) {
        ref++; alt++; pos++;
        reflen--; altlen--;
      }

      // Right trim
      while(reflen && altlen && ref[reflen-1] == alt[altlen-1]) {
        reflen--; altlen--;
      }

      // Only deal with insertions / deletions
      if(reflen != altlen)
      {
        // Fetch window bp either side
        int chromlen = faidx_seq_len(fai, bcf_seqname(hdrin, v));
        // start end coordinates are inclusive
        int start = MAX2(pos-window, 0);
        int end = MIN2(pos+reflen+window, chromlen-1);
        int len = 0;
        char *reg = faidx_fetch_seq(fai, bcf_seqname(hdrin, v), start, end, &len);
        if(!reg) die("Cannot find ref: %s", bcf_seqname(hdrin, v));
        if(len != end-start+1) die("Region not as expected %i vs %i", len, end-start+1);

        // calculate
        hp = calc_homopolymer(ref, reflen, alt, altlen, pos-start, reg, len);
        free(reg);

        // Annotate
        if(hp > 1) {
          int a = bcf_update_info_int32(hdrout, v, "HRun", &hp, 1);
          if(a < 0) die("VCF annotation error");
        }
      }
    }

    // Print
    if(bcf_write(vout, hdrout, v) != 0) die("Cannot write record");
  }

  bcf_destroy(v);
  bcf_hdr_destroy(hdrin);
  bcf_hdr_destroy(hdrout);
  hts_close(vin);
  hts_close(vout);

  fai_destroy(fai);

  return EXIT_SUCCESS;
}
