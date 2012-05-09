#ifndef __ALN_H__
#define __ALN_H__

#include "seq.h"

#define ALN_MAX_SEQ_LEN 1024
#define ALN_MAX_SEQ_LEN_PLUS1 1025
#define ALN_DEFAULT_MATCH_SCORE 1
#define ALN_DEFAULT_MISMATCH_SCORE -3
#define ALN_DEFAULT_OTHER_SCORE 0
#define ALN_DEFAULT_GAP_SCORE -5
#define ALN_DEFAULT_MASK_CHAR 'X'
#define ALN_DEFAULT_SCORE_CUTOFF 5
#define ALN_DEFAULT_MAX_QUAL 40

#define ALN_UNDEF_SCORE (-9999999)
#define ALN_MIN_SCORE (ALN_UNDEF_SCORE+1)

typedef struct AlnNode_t AlnNode;

struct AlnNode_t {
  int i; /* row num */
  int j; /* col num */
  int score;

  /* length of total path */
  int path_len;

  /* where did highest-weight path start? */
  int i_start;
  int j_start;

  /* trace back for highest weight path */
  AlnNode *back_ptr;
};


AlnNode *aln_local(AlnNode **aln_matrix, 
		   int **score_matrix,
		   const int gap_score,
		   Seq *seq1, Seq *seq2);

AlnNode *aln_sgbl_end1_start2(AlnNode **aln_matrix, 
			      int **score_matrix,
			      const int gap_score,
			      Seq *seq1, Seq *seq2);

void aln_matrix_free(AlnNode **matrix);
AlnNode **aln_matrix_new();


int **aln_score_matrix_new(const int match_score, const int mismatch_score, 
			   const int other_score);


void aln_score_matrix_free(int **score_matrix);

void aln_write(FILE *f, AlnNode *end, Seq *seq1, Seq *seq2);


void aln_get_nucs(AlnNode *end, Seq *seq1, Seq *seq2,
		  const unsigned char *qual1,
		  const unsigned char *qual2,
		  unsigned char *nuc_buf1,
		  unsigned char *nuc_buf2,
		  unsigned char *qual_buf1,
		  unsigned char *qual_buf2,
		  const size_t buf_sz);

#endif
