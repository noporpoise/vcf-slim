#ifndef COMMON_H_

#include "hts.h"
#include "vcf.h"
#include "faidx.h"

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MAX2(a,b) ((a) >= (b) ? (a) : (b))

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0)

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
  // Check ref allele is within the reference and matches
  if(v->pos + v->rlen > *chromlen) {
    die("Ref allele goes out of bounds: %s %i %s %s [chromlen: %i]",
        bcf_seqname(hdr, v), v->pos+1, v->d.id, v->d.allele[0], *chromlen);
  }
}

// Trim bases that match with the ref
static inline size_t trimmed_alt(const bcf1_t *v, size_t aid,
                                 size_t *rptr, size_t *aptr)
{
  const char *ref = v->d.allele[0], *alt = v->d.allele[aid];
  size_t rshift = 0;
  size_t reflen = v->rlen;
  size_t altlen = strlen(alt);

  // Left trim
  while(reflen && altlen && *ref == *alt) {
    ref++; alt++; rshift++;
    reflen--; altlen--;
  }

  // Right trim
  while(reflen && altlen && ref[reflen-1] == alt[altlen-1]) {
    reflen--; altlen--;
  }

  *rptr = reflen;
  *aptr = altlen;

  return rshift;
}

static inline void trim_alleles(const bcf1_t *v, size_t *ltrimptr, size_t *rtrimptr)
{
  size_t i, ltrim, rtrim, tmp_ltrim, tmp_rtrim, rlen2, alen2;
  ltrim = trimmed_alt(v, 1, &rlen2, &alen2);
  rtrim = v->rlen - (ltrim + rlen2);

  for(i = 2; i < (size_t)v->n_allele; i++) {
    tmp_ltrim = trimmed_alt(v, i, &rlen2, &alen2);
    tmp_rtrim = v->rlen - (tmp_ltrim + rlen2);
    ltrim = MIN2(ltrim, tmp_ltrim);
    rtrim = MIN2(rtrim, tmp_rtrim);
  }

  // fprintf(stderr, " ltrim:%zu rtrim:%zu\n", ltrim, rtrim);

  *ltrimptr = ltrim;
  *rtrimptr = rtrim;
}


#endif /* COMMON_H_ */
