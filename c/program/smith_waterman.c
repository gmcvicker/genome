#include <stdio.h>
#include <zlib.h>

#include "aln.h"
#include "seq.h"
#include "fastq.h"
#include "util.h"
#include "memutil.h"




int main(int argc, char **argv) {
  gzFile q_fastq_gz, t_fasta_gz;
  FastqSeq read;
  Seq *t_seq, *q_seq;
  int status;
  AlnNode **matrix, *end;
  int **score_matrix;
  long n_row, n_col;
  long n_reads, total_score, n_pos_score, n_neg_score;
  double mean_score;

  if(argc != 3) {
    fprintf(stderr, "usage: %s query_seqs1.fq target_seqs2.fa\n", argv[0]);
    exit(2);
  }

  score_matrix = aln_score_matrix_new(ALN_DEFAULT_MATCH_SCORE,
				      ALN_DEFAULT_MISMATCH_SCORE,
				      ALN_DEFAULT_OTHER_SCORE);

  t_seq = seq_new();
  q_seq = seq_new();
  
  t_fasta_gz = util_must_gzopen(argv[2], "rb");


  n_row = 1000;
  n_col = 1000;
  matrix = aln_matrix_new(n_row, n_col);
  
  
  while((seq_read_fasta_record(t_seq, t_fasta_gz)) >= 0) {

    fprintf(stderr, "--------\n%s\n", t_seq->name);    
    q_fastq_gz = util_must_gzopen(argv[1], "rb");

    /* expand alignment matrix if needed */
    if(n_col < t_seq->len) {
      n_col = t_seq->len * 2;
      fprintf(stderr, "expanding alignment matrix to %ldx%ld\n", n_row, n_col);
      aln_matrix_free(matrix);
      matrix = aln_matrix_new(n_row, n_col);
    }

    n_reads = 0;
    total_score = 0;
    n_pos_score = 0;
    n_neg_score = 0;
    
    while((status = fastq_parse_read(&read, q_fastq_gz)) != FASTQ_END) {
      if(status == FASTQ_ERR) {
	my_warn("%s:%d: error parsing read", __FILE__, __LINE__);
      }
      
      seq_read_seqstr(q_seq, read.line2);
    
      /* expand alignment matrix if needed */
      if(n_row < q_seq->len) {
	n_row = q_seq->len * 2;
	fprintf(stderr, "expanding alignment matrix to %ldx%ld\n",
		n_row, n_col);
	aln_matrix_free(matrix);
	matrix = aln_matrix_new(n_row, n_col);
      }
      
      end = aln_semiglobal(matrix, score_matrix, ALN_DEFAULT_GAP_SCORE,
			   q_seq, t_seq);

      n_reads += 1;
      total_score += end->score;
      if(end->score >= 0) {
	n_pos_score += 1;
      } else {
	n_neg_score += 1;
      }
    }
    
    gzclose(q_fastq_gz);

    mean_score = (double)total_score / (double)n_reads;

    fprintf(stdout, "%s\t%ld\t%ld\t%ld\t%ld\t%.2f\n",
	    t_seq->name, n_reads, total_score, n_pos_score, n_neg_score,
	    mean_score);    
  }

  
  fprintf(stderr, "freeing memory\n");
  aln_matrix_free(matrix);
  aln_score_matrix_free(score_matrix);
  seq_free(t_seq);
  seq_free(q_seq);
  gzclose(t_fasta_gz);
}
