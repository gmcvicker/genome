#include <stdio.h>
#include <string.h>

#include "memutil.h"
#include "nuc.h"
#include "util.h"
#include "err.h"

static const char NUC_SYMBOL[NUM_NUCS] = {'A', 'C', 'G', 'T', '-', 'N'};


/**
 * Converts a nucleotide ID into a character.
 */
char nuc_id_to_char(const unsigned char id) {
  if(id > NUM_NUCS) {
    my_err("Invalid nucleotide identifier %d", id);
  }
  return NUC_SYMBOL[id];
}



/**
 * Converts a nucleotide to a unique nucleotide integer identifier.
 */
unsigned char nuc_char_to_id(const char nuc) {
  switch(nuc) {
  case('A'): case('a'): return NUC_A;
  case('C'): case('c'): return NUC_C;
  case('T'): case('t'): return NUC_T;
  case('G'): case('g'): return NUC_G;
  case('.'): case('-'): case('*'): return NUC_GAP;
  }
  return NUC_N;
}



/**
 * Fills the provided buffer with a null-terminated string representation 
 * of the provided array of nucleotide arrays. The provided buffer
 * must be at least len+1 bytes long, and the provided nucleotide
 * id array must be at least len bytes long.
 *
 * If the provided buffer is NULL, a new buffer of length len+1 is
 * allocated and returned.
 */
char *nuc_ids_to_str(char *buf, const unsigned char *ids, const long len) {
  long i;
  
  if(buf == NULL) {
    buf = my_new(char, len+1);
  }

  for(i = 0; i < len; i++) {
    buf[i] = nuc_id_to_char(ids[i]);
  }
  buf[len] = '\0';
  
  return buf;
}



/**
 * Fills the provided buffer with nucleotide ids corresponding to the
 * string representation of DNA sequence provided. The returned array
 * is NOT null terminated. The len argument should correspond to the
 * length of the provided DNA character string. The provided buf must
 * be at least len bytes long. If the provided buf is NULL a new buffer
 * of length len is allocated and returned.
 *
 * The return value is just the ptr to the provided buf.
 */
unsigned char *nuc_str_to_ids(unsigned char *buf, const char *str, 
			      const long len) {
  long i;
  
  if(buf == NULL) {
    buf = my_new(unsigned char, len);
  }

  for(i = 0; i < len; i++) {
    buf[i] = nuc_char_to_id(str[i]);
  }
  return buf;
}

