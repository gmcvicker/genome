
#ifndef __FASTA_H__
#define __FASTA_H__

#include <zlib.h>

#define FASTA_BUF_SZ 1024

typedef struct {
  char *path;
  char *header;
  char *seqstr;
  long seqlen;
} FASTA;


void fasta_read_record(gzFile *f, FASTA *fasta);
FASTA *fasta_read_file_array(char *filename, long *num_read);
FASTA *fasta_read_file(char *filename);
void fasta_free(FASTA *fasta);

#endif
