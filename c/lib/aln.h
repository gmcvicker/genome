#ifndef __ALN_H__
#define __ALN_H__

#include "seq.h"

#define ALN_DEFAULT_MATCH_SCORE 1
#define ALN_DEFAULT_MISMATCH_SCORE -1
#define ALN_DEFAULT_OTHER_SCORE 0
#define ALN_DEFAULT_GAP_OPEN_SCORE -2
#define ALN_DEFAULT_GAP_EXT_SCORE -1

#define ALN_UNDEF_SCORE (-9999999)


#define ALN_TYPE_GAP1 1 /* gap in sequence 1 */
#define ALN_TYPE_GAP2 2 /* gap in sequence 2 */
#define ALN_TYPE_MM 5 /* match or mismatch */

typedef struct AlnNode_t AlnNode;

struct AlnNode_t {
  long i; /* row num */
  long j; /* col num */
  long n_row; /* total number of rows in matrix */
  long n_col; /* total number of cols in matrix */
  long score;

  int type;

  /* length of total path */
  long path_len;

  /* where did highest-weight path start? */
  long i_start;
  long j_start;

  /* trace back for highest weight path */
  AlnNode *back_ptr;
};


AlnNode *aln_local(AlnNode **aln_matrix,
		   int **score_matrix,
		   int gap_open, int gap_ext,
		   Seq *seq1, Seq *seq2);


AlnNode *aln_semiglobal(AlnNode **aln_matrix,
			int **score_matrix,
			const int gap_score,
			Seq *seq1, Seq *seq2);


AlnNode *aln_semiglobal_end1_start2(AlnNode **aln_matrix,
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
		  unsigned char *qual_buf2);

#endif
