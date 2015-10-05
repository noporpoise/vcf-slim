#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>

#include "hts.h"
#include "vcf.h"
#include "faidx.h"

const bool TRIM_ALLELES = false;
const char *cmd = NULL;

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0)

static void print_usage()
{
  fprintf(stderr, "usage: %s <dist> <in.vcf>\n", cmd);
  fprintf(stderr, "  Filter out entries within <dist> bp of each other.\n");
  fprintf(stderr, "  Input must be sorted. We trim alleles. Print to STDOUT.\n");
  exit(-1);
}

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MAX2(a,b) ((a) >= (b) ? (a) : (b))

#define keepv(a,b,dist) ((a)->rid != (b)->rid || (a)->pos + (a)->rlen + (dist) <= (b)->pos)
#define write_var(vout,hdr,v) do { \
  if(bcf_write(vout, hdr, v) != 0) die("Cannot write record"); \
} while(0)

static inline int var_left_trim(const bcf1_t *v)
{
  int t, j;
  for(t = 0; t < v->rlen; t++) { // loop over ref allele bases
    for(j = 1; j < v->n_allele && v->d.allele[j][t] == v->d.allele[0][t]; j++)
    if(j < v->n_allele) break;
  }
  return t;
}

static inline int var_right_trim(const bcf1_t *v)
{
  // Loop over alleles
  int i, alen, trim = v->rlen, t;
  for(i = 1; i < v->n_allele && trim; i++) {
    alen = strlen(v->d.allele[i]);
    trim = MIN2(trim, alen);
    t = 0;
    while(t < trim && v->d.allele[i][alen-1-t] == v->d.allele[0][v->rlen-1-t]) { t++; }
    trim = t;
  }
  return trim;
}

// Push up start by matching bases
static inline int var_start(const bcf1_t *v)
{
  return v->pos + (TRIM_ALLELES ? var_left_trim(v) : 0);
}

static inline int var_end(const bcf1_t *v)
{
  return v->pos + v->rlen - (TRIM_ALLELES ? var_right_trim(v) : 0);
}

// Returns true if the distance between a, b is great enough to print
static inline bool keepv2(const bcf1_t *a, const bcf1_t *b, int dist, int *end)
{
  if(a->rid != b->rid) { *end = var_end(b); return true; }
  int bstart = var_start(b);
  int bend = var_end(b);
  int last_end = *end;
  *end = MAX2(*end, bend);
  return last_end + dist <= bstart;
}

int main(int argc, char **argv)
{
  cmd = argv[0];
  if(argc != 3) print_usage();

  int dist = atoi(argv[1]);
  const char *vcf_path = argv[2];
  const char *out_path = "-";

  // open input / output
  htsFile *vin = hts_open(vcf_path, "r");
  if(vin == NULL) die("Cannot read %s: %s", vcf_path, strerror(errno));
  htsFile *vout = hts_open(out_path, "w");
  if(vout == NULL) die("Cannot write %s: %s", out_path, strerror(errno));

  // read/write headers
  bcf_hdr_t *hdr = bcf_hdr_read(vin);
  if(bcf_hdr_write(vout, hdr) != 0) die("Cannot write header");

  // read vcf entries
  size_t i, nv;
  int maxend = 0;

  bcf1_t *v[3];
  v[0] = bcf_init();
  v[1] = bcf_init();
  v[2] = bcf_init();

  for(nv = 0; nv < 3; nv++) {
    if(bcf_read(vin, hdr, v[nv]) < 0) break;
    bcf_unpack(v[nv], BCF_UN_STR);
  }

  if(nv == 1) {
    write_var(vout, hdr, v[0]);
  }
  else if(nv == 2) {
    maxend = var_end(v[0]);
    if(keepv2(v[0], v[1], dist, &maxend)) {
      write_var(vout, hdr, v[0]);
      write_var(vout, hdr, v[1]);
    }
  }
  else if(nv == 3)
  {
    maxend = var_end(v[0]);
    bool kn, kp = keepv2(v[0], v[1], dist, &maxend);
    if(kp) write_var(vout, hdr, v[0]);

    do
    {
      bcf_unpack(v[2], BCF_UN_STR);
      kn = keepv2(v[1], v[2], dist, &maxend);
      if(kp && kn) write_var(vout, hdr, v[1]);
      kp = kn;
      bcf1_t *tmp = v[0]; v[0] = v[1]; v[1] = v[2]; v[2] = tmp;
    } while(bcf_read(vin, hdr, v[2]) >= 0);

    if(kp) write_var(vout, hdr, v[1]);
  }

  for(i = 0; i < 3; i++) bcf_destroy(v[i]);

  bcf_hdr_destroy(hdr);
  hts_close(vin);
  hts_close(vout);

  return EXIT_SUCCESS;
}
