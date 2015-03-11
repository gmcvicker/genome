#include <stdio.h>
#include <stdlib.h>


#include "seq.h"


int main(int argc, char **argv) {
  char *filename;
  long seq_len;
  Seq *seq;
  
  if(argc != 2) {
    fprintf(stderr, "usage: test_fasta <file.fa>\n");
    return -1;
  }

  filename = argv[1];

  fprintf(stderr, "reading from file '%s'\n", filename);

  seq = seq_new();
  seq_len = seq_read_fasta_from_file(seq, filename);

  fprintf(stderr, "read fasta record with %ld bp of sequence\n", seq_len);
  
  seq_free(seq);

  fprintf(stderr, "done\n");
  
  return 0;
}
