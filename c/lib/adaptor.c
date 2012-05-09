
#include <zlib.h>
#include <string.h>

#include "adaptor.h"
#include "seq.h"
#include "err.h"
#include "util.h"

Seq *adaptor_read_seq(const char *filename) {
  Seq *adp_seq;
  gzFile *adp_file;

  if(filename == NULL) {
    return NULL;
  }

  /* read adaptor sequence */
  fprintf(stderr, "adaptor file=%s\n", filename);

  adp_seq = seq_new();
  adp_file = gzopen(filename, "r");
  if(!adp_file) {
    my_err("%s:%d: could not open adaptor file: %s", __FILE__, __LINE__,
	  filename);
  }
  if(!seq_read_fasta_record(adp_seq, adp_file)) {
    my_err("%s:%d: could not read fasta record from adaptor file: %s",
	  __FILE__, __LINE__, filename);
  }

  return adp_seq;
}


/**
 * Clears the mask for the provided sequence
 */
void adaptor_mask_clear(Seq *seq, unsigned char *mask) {
  memset(mask, 0, seq->len);
}


/**
 * masks the left end of the sequence given the provided
 * alignment to the left adaptor
 */
void adaptor_mask_left(Seq *seq, unsigned char *mask,
		       AlnNode *aln_end) {
  
  long len;

  if(aln_end->j_start != 0) {
    my_err("%s:%d: expected left adaptor alignment to start at 1, not %ld",
	  __FILE__, __LINE__, aln_end->j_start+1);
  }

  len = aln_end->j + 1;
  if(len >= seq->len) {
    my_err("%s:%d: mask coordinates (1-%ld) fall outside of sequence "
	  "range (1-%ld)", __FILE__, __LINE__, len, seq->len);
  }

  memset(mask, TRUE, len);
}



/**
 * masks the right end of the sequence given the provided
 * alignment to the right adaptor
 */
void adaptor_mask_right(Seq *seq, unsigned char *mask,
			AlnNode *aln_end) {
  long start_idx, len;
  
  if(aln_end->i+1 != seq->len) {
    my_err("%s:%d: expected right adaptor alignment to end at %ld, not %ld",
	  __FILE__, __LINE__, seq->len, aln_end->i+1);
  }

  start_idx = aln_end->i_start;
  len = aln_end->i - start_idx + 1;
  memset(&mask[start_idx], TRUE, len);

  if(len + start_idx > seq->len) {
    my_err("%s:%d: mask coordinates (%ld-%ld) fall outside of sequence "
	  "range (1-%ld)", __FILE__, __LINE__, start_idx+1, start_idx+len,
	  seq->len);
  }
}


