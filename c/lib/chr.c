
#include <limits.h>

#include "memutil.h"
#include "chr.h"
#include "util.h"



/**
 * Returns a (deep) copy of the provided chromosome
 */
Chromosome *chr_copy(const Chromosome *chr) {
  Chromosome *new_chr;

  new_chr = my_new(Chromosome, 1);
  new_chr->id = chr->id;

  if(chr->name) {
    new_chr->name = util_str_dup(chr->name);
  } else {
    new_chr->name = NULL;
  }

  if(chr->assembly) {
    new_chr->assembly = util_str_dup(chr->assembly);
  } else {
    new_chr->assembly = NULL;
  }

  new_chr->len = chr->len;
  
  return new_chr;
}


/**
 * Reads an array chromosomes from a file containing a name and length on
 * each line, separated by a white space character.
 */
Chromosome *chr_read_file(const char *filename, int *n_chr) {
  char buf[LINE_MAX], name_buf[LINE_MAX];
  Chromosome *chrs;
  FILE *f;
  int n, i;

  if(util_has_gz_ext(filename)) {
    my_err("%s:%d: gzipped chr files not currently supported\n",
	   __FILE__, __LINE__);
  }

  f = util_must_fopen(filename, "r");
  *n_chr = util_fcount_lines(f);
  
  if(*n_chr < 1) {
    my_err("%s:%d: chromosome file '%s' is empty\n", 
	   __FILE__, __LINE__, filename);
  }

  chrs = my_new(Chromosome, *n_chr);
  for(i = 0; i < *n_chr; i++) {
    if(!fgets(buf, sizeof(buf), f)) {
      my_err("%s:%d: expected %d lines in file, but only read %d\n", 
	     __FILE__, __LINE__, *n_chr, i);
    }

    n = sscanf(buf, "%s %ld", name_buf, &chrs[i].len);
    if(n < 2) {
      my_err("%s:%d: line did not have at least 2 tokens:\n'%s'",
	     __FILE__, __LINE__, buf);
    }
    chrs[i].name = util_str_dup(name_buf);
    chrs[i].assembly = NULL;
    chrs[i].id = i;
    
    if(chrs[i].len < 1) {
      my_err("%s:%d: chr length (%ld) should be >= 1",
	     __FILE__, __LINE__, chrs[i].len);
    }
  }

  fclose(f);

  return chrs;
}



/**
 * Frees memory allocated for an array of chromosomes
 */
void chr_array_free(Chromosome *chrs, int n_chr) {
  int i;

  for(i = 0; i < n_chr; i++) {
    if(chrs[i].name) {
      my_free(chrs[i].name);
    }
    if(chrs[i].assembly) {
      my_free(chrs[i].assembly);
    }
  }

  my_free(chrs);
}


/**
 * Frees memory that was allocated for a single chromosome
 */
void chr_free(Chromosome *chr) {
  if(chr->name != NULL) {
    my_free(chr->name);
  }

  if(chr->assembly != NULL) {
    my_free(chr->assembly);
  }

  my_free(chr);
}

