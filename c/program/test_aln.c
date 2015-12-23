#include <stdio.h>
#include <zlib.h>

#include "aln.h"
#include "seq.h"
#include "fastq.h"
#include "util.h"
#include "memutil.h"


int main(int argc, char **argv) {
  Seq *t_seq, *q_seq;
  int **score_matrix;
  long gap_open, gap_ext, match_score, mismatch_score;
  AlnNode **matrix, *end;

  if(argc != 7) {
    fprintf(stderr, "usage: %s <query_seq> <target_seq> <match_score> <mismatch_score> <gap_open> <gap_ext>\n",
	    argv[0]);
    exit(2);
  }

  t_seq = seq_new();
  q_seq = seq_new();
  seq_read_seqstr(t_seq, argv[1]);
  seq_read_seqstr(q_seq, argv[2]);

  match_score = util_parse_long(argv[3]);
  mismatch_score = util_parse_long(argv[4]);
  gap_open = util_parse_long(argv[5]);
  gap_ext = util_parse_long(argv[6]);

  fprintf(stderr, "match: %ld\n"
	  "mismatch: %ld\n"
	  "gap_open: %ld\n"
	  "gap_ext: %ld\n", match_score, mismatch_score, gap_open, gap_ext);

  score_matrix = aln_score_matrix_new(match_score,
				      mismatch_score,
				      ALN_DEFAULT_OTHER_SCORE);

  matrix = aln_matrix_new(q_seq->len, t_seq->len);

  end = aln_local(matrix, score_matrix, gap_open, gap_ext,
		  q_seq, t_seq);

  aln_write(stdout, end, q_seq, t_seq);

  seq_free(t_seq);
  seq_free(q_seq);
  aln_matrix_free(matrix);
  aln_score_matrix_free(score_matrix);

  return 0;
}
