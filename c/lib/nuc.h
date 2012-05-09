#ifndef __NUC_H__
#define __NUC_H__

#include <stdio.h>
#include <ctype.h>
#include "seq.h"

enum nucleotide {NUC_A=0, NUC_C, NUC_G, NUC_T, NUC_GAP, NUC_N, NUM_NUCS};

/* the number of "true" nucleotides, i.e. not gaps or ambiguity chars */
#define NUM_REAL_NUCS 4


/* macro function, compliments a nucleotide ID */
#define nuc_comp(x) ((x==NUC_A) ? NUC_T : \
                    ((x==NUC_T) ? NUC_A : \
                    ((x==NUC_G) ? NUC_C : \
                    ((x==NUC_C) ? NUC_G : x))))

char nuc_id_to_char(const unsigned char id);
unsigned char nuc_char_to_id(const char nuc);

char *nuc_ids_to_str(char *buf, const unsigned char *ids, const long len);
unsigned char *nuc_str_to_ids(unsigned char *buf, const char *str, 
			      const long len);

double **nuc_read_matrix(const char *filename);

void nuc_write_matrix_dbl(FILE *fh, double **matrix);
void nuc_write_matrix_ul(FILE *fh, long **matrix);

long **nuc_matrix_ul_new(void);
double **nuc_matrix_dbl_new(void);

void nuc_matrix_ul_set_zero(long **matrix);


#endif
