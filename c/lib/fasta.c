#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "err.h"
#include "memutil.h"
#include "fasta.h"
#include "util.h"


/*
 * Helper function, reads a single FASTA record from a file and
 * advances the filehandle to the end of the record.
 */
void fasta_read_record(gzFile *f, FASTA *fasta) {
  long len, seq_pos, i;
  z_off_t fseq_start, fseq_end;
  char buf[FASTA_BUF_SZ];

  /* read the header */
  fasta->header = util_gzgets_line(f);
  fasta->path = NULL;

  if(fasta->header[0] != '>') {
    my_err("fasta_read_record: expected header to begin with '>'");
  }

  /* record current position in file */
  fseq_start = gztell(f);

  /* read ahead in the file to determine length of current sequence */
  fasta->seqlen = 0;
  while(gzgets(f, buf, FASTA_BUF_SZ) != NULL && buf[0] != '>') {
    len = strlen(buf);

    if(len == 0) {
      continue;
    }

    if(buf[len-1] == '\n') {
      fasta->seqlen += len - 1;
    } else {
      fasta->seqlen += len;
    }
  }

  /* allocate room for sequence */
  fasta->seqstr = my_new(char, fasta->seqlen);

  /* Read the sequence again, this time writing to the allocated mem */
  gzseek(f, fseq_start, SEEK_SET);
  fseq_end = gztell(f);
  seq_pos = 0;
  while(gzgets(f, buf, FASTA_BUF_SZ) != NULL) {
    if(buf[0] == '>') {
      /* go back to where we ended last sequence */
      gzseek(f, fseq_end, SEEK_SET);
      break;
    }

    len = strlen(buf);

    if(len == 0) {
      continue;
    }

    if(buf[len-1] == '\n') {
      len -= 1;
      if(len == 0) {
	continue;
      }
    }

    for(i = 0; i < len; i++) {
      fasta->seqstr[seq_pos + i] = buf[i];
    }
    seq_pos += len;
    
    /* remember location, in case next line is new sequence header */
    fseq_end = gztell(f);
  }

  return;
}




/*
 * Reads all of the sequence records from a FASTA file into an
 * array. The number of sequences read is assigned to the value
 * pointed to by the num_read argument. The array of sequences should
 * be freed when no longer needed.
 */
FASTA *fasta_read_file_array(char *filename, long *num_read) {
  FASTA *records;
  long i;
  gzFile *f;

  f = gzopen(filename, "rb");
  if(f == NULL) {
    my_err("fasta_read_file_array: could not open file '%s'", filename);
  }


  /* count the number of FASTA header lines in the file */
  *num_read = util_gzcount_lines_match(f, ">");

  records = my_new(FASTA, *num_read);
  for(i = 0; i < *num_read; i++) {
    fasta_read_record(f, &records[i]);
    records[i].path = util_str_dup(filename);
  }

  gzclose(f);
  
  return records;
}



/**
 * Reads a single sequence from a FASTA file.  Terminates with an
 * error if the number of sequences in the file is not 1.
 */
FASTA *fasta_read_file(char *filename) {
  FASTA *record;
  long num;

  record = fasta_read_file_array(filename, &num);

  if(num != 1) {
    my_err("%s:%d: expected file '%s' to contain 1 FASTA record"
	    ", but contains %ld.\nThe first line should start with '>'.", 
	    __FILE__, __LINE__, filename, num);
  }

  return record;
}




/**
 * Frees memory allocated for a FASTA data structure
 */
void fasta_free(FASTA *fasta) {
  if(fasta->seqstr != NULL) {
    my_free(fasta->seqstr);
  }

  if(fasta->header != NULL) {
    my_free(fasta->header);
  }

  if(fasta->path != NULL) {
    my_free(fasta->path);
  }

  my_free(fasta);
}
