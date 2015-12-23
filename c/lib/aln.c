
#include "err.h"
#include "aln.h"
#include "nuc.h"
#include "memutil.h"
#include "util.h"
#include "seq.h"


int **aln_score_matrix_new(const int match_score, const int mismatch_score,
			   const int other_score) {
  unsigned char i, j;
  int **score_matrix;

  score_matrix = my_new(int *, NUM_NUCS);

  for(i = 0; i < NUM_NUCS; i++) {
    score_matrix[i] = my_new(int, NUM_NUCS);

    for(j = 0; j < NUM_NUCS; j++) {
      if(i == NUC_N || i == NUC_GAP ||
         j == NUC_N || j == NUC_GAP) {
	score_matrix[i][j] = other_score;
      }
      else if(i == j) {
	score_matrix[i][i] = match_score;
      }
      else {
	score_matrix[i][j] = mismatch_score;
      }
    }
  }

  return score_matrix;
}



void aln_score_matrix_free(int **score_matrix) {
  unsigned char i;

  for(i = 0; i < NUM_NUCS; i++) {
    my_free(score_matrix[i]);
  }
  my_free(score_matrix);
}




/**
 * Creates new alignment node matrix with dimension n_row x n_col
 */
AlnNode **aln_matrix_new(long n_row, long n_col) {
  AlnNode **matrix;
  long i, j;

  if(n_row < 1) {
    my_err("%s:%d: n_row must be at least 1",
	   __FILE__, __LINE__);
  }
  if(n_col < 1) {
    my_err("%s:%d: n_col must be at least 1",
	   __FILE__, __LINE__);
  }

  matrix = my_new(AlnNode *, n_row);

  for(i = 0; i < n_row; i++) {
    matrix[i] = my_new(AlnNode, n_col);
    for(j = 0; j < n_col; j++) {
      matrix[i][j].n_row = n_row;
      matrix[i][j].n_col = n_col;
      matrix[i][j].i = i;
      matrix[i][j].j = j;
      matrix[i][j].score = ALN_UNDEF_SCORE;
      matrix[i][j].path_len = 0;
      matrix[i][j].back_ptr = NULL;
    }
  }

  return matrix;
}



void aln_matrix_free(AlnNode **matrix) {
  long i, n_row;

  n_row = matrix[0][0].n_row;

  for(i = 0; i < n_row; i++) {
    my_free(matrix[i]);
  }
  my_free(matrix);
}



/**
 * Performs a local aligment of the two provided sequences
 */
AlnNode *aln_local(AlnNode **aln_matrix, int **score_matrix,
		   int gap_open, int gap_ext, Seq *seq1, Seq *seq2) {
  long i, j;
  long n_row;
  long n_col;
  int mm_score, gap_score;
  AlnNode *prev, *cur, *max_node;

  n_row = aln_matrix[0][0].n_row;
  n_col = aln_matrix[0][0].n_col;

  /* TODO: here we could expand the alignment matrix as necessary
   * rather than giving an error
   */
  if(seq1->len > n_row) {
    my_err("%s:%d: seq1 sequence length exceeds number of rows (%ld)",
	    __FILE__, __LINE__, n_row);
  }
  if(seq2->len > n_col) {
    my_err("%s:%d: seq2 sequence length exceeds numer of cols (%ld)",
	    __FILE__, __LINE__, n_col);
  }

  max_node = NULL;

  for(i = 0; i < seq1->len; i++) {
    for(j = 0; j < seq2->len; j++) {
      /* get match/mismatch score for this position */
      cur = &aln_matrix[i][j];
      mm_score = score_matrix[seq1->sym[i]][seq2->sym[j]];

      /* Score for starting new local alignment at this position
       * is just match/mistmatch score
       */
      cur->score    = mm_score;
      cur->type     = ALN_TYPE_MM;
      cur->back_ptr = NULL;
      cur->i_start  = i;
      cur->j_start  = j;
      cur->path_len = 1;

      /* does inserting GAP from ABOVE (in seq1) give higher score? */
      if(i > 0) {
	prev = &aln_matrix[i-1][j];

	gap_score = (prev->type == ALN_TYPE_GAP1) ? gap_ext : gap_open;

	if((prev->score + gap_score) > cur->score) {
	  /* extend existing gap */
	  cur->score    = prev->score + gap_score;
	  cur->back_ptr = prev;
	  cur->i_start  = prev->i_start;
	  cur->j_start  = prev->j_start;
	  cur->path_len = prev->path_len + 1;
	  cur->type = ALN_TYPE_GAP1;
	}
      }

      /* does inserting GAP from LEFT (gap in seq2) give higher score? */
      if(j > 0) {
	prev = &aln_matrix[i][j-1];

	gap_score = (prev->type == ALN_TYPE_GAP2) ? gap_ext : gap_open;

	if((prev->score + gap_score) > cur->score) {
	    /* extend existing gap */
	    cur->score    = prev->score + gap_score;
	    cur->back_ptr = prev;
	    cur->i_start  = prev->i_start;
	    cur->j_start  = prev->j_start;
	    cur->path_len = prev->path_len + 1;
	    cur->type     = ALN_TYPE_GAP2;
	}
      }

      /* does extending from ABOVE,LEFT gives higher score? */
      if(i > 0 && j > 0) {
	prev = &aln_matrix[i-1][j-1];
	if((prev->score + mm_score) > cur->score) {
	  cur->score    = prev->score + mm_score;
	  cur->back_ptr = prev;
	  cur->i_start  = prev->i_start;
	  cur->j_start  = prev->j_start;
	  cur->path_len = prev->path_len + 1;
	  cur->type = ALN_TYPE_MM;
	}
      }

      if((max_node == NULL) || (cur->score > max_node->score)) {
	max_node = cur;
      }
    }
  }

  return max_node;
}


/**
 * Performs semi-global alignment of sequence 1 against
 * sequence 2. I.e. the full length of sequence1 must be
 * spanned by the alignment, but not of sequence2.
 */
AlnNode *aln_semiglobal(AlnNode **aln_matrix, int **score_matrix,
			 const int gap_score, Seq *seq1, Seq *seq2) {
  long i, j;
  long n_row;
  long n_col;
  int mm_score;
  AlnNode *prev, *cur, *max_node;

  /** TODO: separate GAP open / ext **/

  n_row = aln_matrix[0][0].n_row;
  n_col = aln_matrix[0][0].n_col;

  /* TODO: here we could expand the alignment matrix as necessary
   * rather than giving an error
   */
  if(seq1->len > n_row) {
    my_err("%s:%d: seq1 sequence length exceeds number of rows (%ld)",
	    __FILE__, __LINE__, n_row);
  }
  if(seq2->len > n_col) {
    my_err("%s:%d: seq2 sequence length exceeds numer of cols (%ld)",
	    __FILE__, __LINE__, n_col);
  }

  /* initialize all columns for first row, this is at beginning of seq1
   * and is where alignment must start
   */
  for(j = 0; j < seq2->len; j++) {
    /* get match/mismatch score for this position */
      cur = &aln_matrix[0][j];
      mm_score = score_matrix[seq1->sym[0]][seq2->sym[j]];

      cur->score    = mm_score;
      cur->back_ptr = NULL;
      cur->i_start  = 0;
      cur->j_start  = j;
      cur->path_len = 1;
      cur->type = ALN_TYPE_MM;
  }

  for(i = 1; i < seq1->len; i++) {
    for(j = 0; j < seq2->len; j++) {
      /* get match/mismatch score for this position */
      cur = &aln_matrix[i][j];
      mm_score = score_matrix[seq1->sym[i]][seq2->sym[j]];

      /* three possible options:
       *   insert GAP from ABOVE,
       *   insert GAP from LEFT,
       *   extending from ABOVE LEFT
       */

      /* GAP from ABOVE */
      prev = &aln_matrix[i-1][j];
      cur->score = prev->score + gap_score;

      cur->back_ptr = prev;
      cur->i_start  = prev->i_start;
      cur->j_start  = prev->j_start;
      cur->path_len = prev->path_len + 1;
      cur->type = ALN_TYPE_GAP1;


      /* does inserting GAP from LEFT give better score? */
      if(j > 0) {
	prev = &aln_matrix[i][j-1];
	if((prev->score + gap_score) > cur->score) {
	  cur->score    = prev->score + gap_score;
	  cur->back_ptr = prev;
	  cur->i_start  = prev->i_start;
	  cur->j_start  = prev->j_start;
	  cur->path_len = prev->path_len + 1;
	  cur->type = ALN_TYPE_GAP2;
	}

	/* does extending from ABOVE, LEFT gives higher score? */
	prev = &aln_matrix[i-1][j-1];
	if((prev->score + mm_score) > cur->score) {
	  cur->score    = prev->score + mm_score;
	  cur->back_ptr = prev;
	  cur->i_start  = prev->i_start;
	  cur->j_start  = prev->j_start;
	  cur->path_len = prev->path_len + 1;
	  cur->type = ALN_TYPE_MM;
	}
      }
    }
  }

  /* we only set max node in last column, since alignment must
   * extend to end of sequence 1
   */
  i = seq1->len-1;
  max_node = &aln_matrix[i][0];
  for(j = 1; j < seq2->len; j++) {
    cur = &aln_matrix[i][j];
    if(cur->score > max_node->score) {
      max_node = cur;
    }
  }

  return max_node;
}



/**
 * Performs semi-global alignment, constrained to the
 * end of sequence 1 and the beginning of sequence 2.
 */
AlnNode *aln_semiglobal_end1_start2(AlnNode **aln_matrix,
				    int **score_matrix,
				    const int gap_score,
				    Seq *seq1, Seq *seq2) {
  long i, j, n_row, n_col;
  int new_score, max_score, mm_score;
  AlnNode *node, *max_node;

  /* TODO: here we could expand the alignment matrix as necessary
   * rather than giving an my_err
   */
  n_row = aln_matrix[0][0].n_row;
  n_col = aln_matrix[0][0].n_col;
  if(seq1->len > n_row) {
    my_err("%s:%d: seq1 sequence length exceeds max sequence length (%ld)",
	    __FILE__, __LINE__, n_row);
  }
  if(seq2->len > n_col) {
    my_err("%s:%d: seq2 sequence length exceeds max sequence length (%ld)",
	    __FILE__, __LINE__, n_col);
  }

  /* number read bases i, adaptor bases j */

  /* we must start at beginning of adaptor sequence where j = 0 */
  /* loop over first column (j=0), assigning scores */
  j = 0;
  for(i = 0; i < seq1->len; i++) {
    aln_matrix[i][j].back_ptr = NULL;
    mm_score = score_matrix[seq1->sym[i]][seq2->sym[j]];
    aln_matrix[i][j].score = mm_score;

    aln_matrix[i][j].path_len = 1;
    aln_matrix[i][j].i_start = i;
    aln_matrix[i][j].j_start = j;
  }

  /* currently there is just one complete alignment, so it
   * has the highest score (first adaptor base aligned
   * to last read base)
   */
  max_node = &aln_matrix[seq1->len-1][0];
  max_score = max_node->score;

  /* starting at second column (j=1), fill in score matrix */
  for(i = 0; i < seq1->len; i++) {
    for(j = 1; j < seq2->len; j++) {

      node = &aln_matrix[i][j];
      node->score = ALN_UNDEF_SCORE;

      if(i > 0) {

	/* get score obtained by adding match/mismatch from i-1,j-1 */
	mm_score = score_matrix[seq1->sym[i]][seq2->sym[j]];
	new_score = aln_matrix[i-1][j-1].score + mm_score;
	if(new_score > node->score) {
	  node->score = new_score;
	  node->back_ptr = &aln_matrix[i-1][j-1];

	  /* extend path */
	  node->path_len = aln_matrix[i-1][j-1].path_len + 1;
	  node->i_start = aln_matrix[i-1][j-1].i_start;
	  node->j_start = aln_matrix[i-1][j-1].j_start;
	}

	/* get score obtained by adding gap from i-1,j */
	new_score = aln_matrix[i-1][j].score + gap_score;
	if(new_score > node->score) {
	  node->score = new_score;
	  node->back_ptr = &aln_matrix[i-1][j];

	  node->path_len = aln_matrix[i-1][j].path_len + 1;
	  node->i_start = aln_matrix[i-1][j].i_start;
	  node->j_start = aln_matrix[i-1][j].j_start;
	}
      }

      /* get score obtained from adding gap from i,j-1 */
      new_score = aln_matrix[i][j-1].score + gap_score;
      if(new_score > node->score) {
	node->score = new_score;
	node->back_ptr = &aln_matrix[i][j-1];

	node->path_len = aln_matrix[i][j-1].path_len + 1;
	node->i_start = aln_matrix[i][j-1].i_start;
	node->j_start = aln_matrix[i][j-1].j_start;
      }

      if(i == (seq1->len-1)) {
	/* Last base of read, where we force alignment to end.
	 * Check if this is max scoring path.
	 */
	if(node->score > max_score) {
	  max_score = node->score;
	  max_node = node;
	}
      }
    }
  }

  return max_node;
}



void aln_get_nucs(AlnNode *end,
		  Seq *seq1, Seq *seq2,
		  const unsigned char *qual1,
		  const unsigned char *qual2,
		  unsigned char *nuc_buf1,
		  unsigned char *nuc_buf2,
		  unsigned char *qual_buf1,
		  unsigned char *qual_buf2) {
  AlnNode *cur, *next;
  int idx, use_qual;

  /* use quality scores if quality args are not NULL */
  use_qual = qual1 && qual2 && qual_buf1 && qual_buf2;

  idx = end->path_len-1;
  next = end;
  cur = next->back_ptr;


  /* move backwards through alignment filling in nucleotides */
  while(cur != NULL) {
    if(idx < 0) {
      my_err("%s:%d: alignment is longer than expected", __FILE__, __LINE__);
    }

    if(cur->i < next->i) {
      nuc_buf1[idx] = seq1->sym[next->i];
    } else {
      nuc_buf1[idx] = NUC_GAP;
    }

    if(cur->j < next->j) {
      nuc_buf2[idx] = seq2->sym[next->j];
    } else {
      nuc_buf2[idx] = NUC_GAP;
    }

    if(use_qual) {
      /* also track quality scores */
      if(cur->i < next->i) {
	qual_buf1[idx] = qual1[next->i];
      } else {
	qual_buf1[idx] = 0;
      }

      if(cur->j < next->j) {
	qual_buf2[idx] = qual2[next->j];
      } else {
	qual_buf2[idx] = 0;
      }
    }

    next = cur;
    cur = cur->back_ptr;
    idx--;
  }

  if(idx > 0) {
    my_err("%s:%d: alignment is shorter than expected", __FILE__, __LINE__);
  }
  fprintf(stderr, "idx: %d\n", idx);

  /* fill in last pair of nucleotides as match or mismatch */
  if(next && idx >= 0) {
    nuc_buf1[idx] = seq1->sym[next->i];
    nuc_buf2[idx] = seq2->sym[next->j];
    if(use_qual) {
      qual_buf1[idx] = qual1[next->i];
      qual_buf2[idx] = qual2[next->j];
    }
  }

}


void aln_write(FILE *f, AlnNode *end, Seq *seq1, Seq *seq2) {
  char *nuc_str1, *nuc_str2, *match_str;
  unsigned char *aln_nuc1;
  unsigned char *aln_nuc2;
  int i;

  nuc_str1 = my_new(char, end->path_len+1);
  nuc_str2 = my_new(char, end->path_len+1);
  match_str = my_new(char, end->path_len+1);

  nuc_str1[end->path_len] = '\0';
  nuc_str2[end->path_len] = '\0';
  match_str[end->path_len] = '\0';

  /** TODO: for efficiency could change this function
   * to take buffers, rather than re-alloc'ing mem each time
   */
  aln_nuc1 = my_malloc(end->path_len);
  aln_nuc2 = my_malloc(end->path_len);

  aln_get_nucs(end, seq1, seq2, NULL, NULL,
	       aln_nuc1, aln_nuc2, NULL, NULL);

  for(i = 0; i < end->path_len; i++) {
    if(aln_nuc1[i] == aln_nuc2[i] &&
       aln_nuc1[i] != NUC_GAP &&
       aln_nuc1[i] != NUC_N) {
      match_str[i] = '|';
    } else {
      match_str[i] = ' ';
    }
  }

  nuc_ids_to_str(nuc_str1, aln_nuc1, end->path_len);
  nuc_ids_to_str(nuc_str2, aln_nuc2, end->path_len);

  fprintf(f, "ALIGN: %s score=%ld len=%ld\n", seq1->name,
	  end->score, end->path_len);
  fprintf(f, "%s\n%s\n%s\n\n", nuc_str1, match_str, nuc_str2);

  my_free(aln_nuc1);
  my_free(aln_nuc2);
  my_free(nuc_str1);
  my_free(nuc_str2);
}

