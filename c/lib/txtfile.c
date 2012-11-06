
#include <limits.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>

#include "txtfile.h"
#include "memutil.h"
#include "err.h"


/**
 * splits a line into tokens, returns non-zero on error
 */
static int split_line(char *line, char **tokens, int n_tok) {  
  int i = 0;
  while((tokens[i] = strsep(&line, " \t")) != NULL) {
    if(tokens[i][0] != '\0') {
      i++;
      if(i >= n_tok) {
	break;
      }
    }
  }
  
  if(i != n_tok) {
    my_err("%s:%d: expected at least %d tokens but got %d.",
	  __FILE__, __LINE__, n_tok, i);
    return TXTFILE_ERR;
  }

  return 0;
}


/**
 * Reads int8 values for an entire chromosome from a tab- or
 * space-delimited text file. The column indices of the position and
 * value columns in the text file must be provided. Unspecified values
 * are set to 0. Returns NULL on error.
 */
char *txtfile_read_int8(const char *filename, const long chr_len, 
			const int pos_idx, const int val_idx) {
  char line[TXTFILE_MAX_LINE];
  gzFile gzf;
  char *vals;
  long pos, val, count;
  int n_tok;
  char **tokens;

  if(pos_idx < 0 || val_idx < 0) {
    my_warn("%s:%d: invalid position (%d) and/or value (%d) indices", 
	    __FILE__, __LINE__, pos_idx, val_idx);
    return NULL;
  }
  if(chr_len < 1) {
    my_warn("%s:%d: invalid chromosmoe len (%ld)", 
	    __FILE__, __LINE__, chr_len);
    return NULL;
  }

  /* allocate space to hold first n line tokens, where n is smallest
   * of position or value index
   */
  n_tok = (pos_idx > val_idx) ? pos_idx + 1 : val_idx + 1;
  tokens = my_malloc(sizeof(char *) * n_tok);
  
  /* allocate bytes and set to 0 */
  vals = my_malloc(sizeof(char) * chr_len);
  bzero(vals, chr_len);
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
    my_free(tokens);
    my_free(vals);
    return NULL;
  }

  pos = 1;
  count = 0;
  while(gzgets(gzf, line, sizeof(line))) {
    if(split_line(line, tokens, n_tok)) {
      /* an error occured splitting the line */
      my_free(tokens);
      my_free(vals);
      return NULL;
    }
  
    errno = 0;
    pos = strtol(tokens[pos_idx], NULL, 10);
    if(pos == 0 && errno) {
      my_warn("%s:%d: could not parse position in col %d as long int: '%s'\n", 
	    __FILE__, __LINE__, pos_idx, tokens[pos_idx]);
      my_free(tokens);
      my_free(vals);
      return NULL;
    }

    errno = 0;
    val = strtol(tokens[val_idx], NULL, 10);
    if(val == 0 && errno) {
      my_warn("%s:%d: could not parse value in col %d as long int: '%s'\n",
	      __FILE__, __LINE__, val_idx, tokens[val_idx]);
      my_free(tokens);
      my_free(vals);
      return NULL;
    }

    if(val > CHAR_MAX) {
      my_warn("%s:%d: setting value at position %ld (%ld) to maximum "
	      "allowed for int8 (%d)\n", __FILE__, __LINE__, 
	      pos, val, CHAR_MAX);
      val = CHAR_MAX;
    }

    if((pos > chr_len) || (pos < 1)) {
      my_warn("%s:%d: skipping pos %ld "
	      "(outside chromosome range 1-%ld)\n", __FILE__, __LINE__,
	      pos, chr_len);
    } else {
      vals[pos-1] = (char)val;
      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }
    }
  }

  my_free(tokens);
  fprintf(stderr, "\n");
  gzclose(gzf);
  
  return vals;
}






/**
 * Reads int16 values for an entire chromosome from a tab- or
 * space-delimited text file. The column indices of the position and
 * value columns in the text file must be provided. Unspecified values
 * are set to 0. Returns NULL on error.
 */
short *txtfile_read_int16(const char *filename, const long chr_len, 
			const int pos_idx, const int val_idx) {
  char line[TXTFILE_MAX_LINE];
  gzFile gzf;
  short *vals;
  long pos, val, count;
  int n_tok;
  char **tokens;

  if(pos_idx < 0 || val_idx < 0) {
    my_warn("%s:%d: invalid position (%d) and/or value (%d) indices", 
	    __FILE__, __LINE__, pos_idx, val_idx);
    return NULL;
  }
  if(chr_len < 1) {
    my_warn("%s:%d: invalid chromosmoe len (%ld)", 
	    __FILE__, __LINE__, chr_len);
    return NULL;
  }

  /* allocate space to hold first n line tokens, where n is smallest
   * of position or value index
   */
  n_tok = (pos_idx > val_idx) ? pos_idx + 1 : val_idx + 1;
  tokens = my_malloc(sizeof(char *) * n_tok);
  
  /* allocate integers and set to 0 */
  vals = my_malloc0(sizeof(short) * chr_len);
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
    my_free(tokens);
    my_free(vals);
    return NULL;
  }

  pos = 1;
  count = 0;
  while(gzgets(gzf, line, sizeof(line))) {
    if(split_line(line, tokens, n_tok)) {
      /* an error occured splitting the line */
      my_free(tokens);
      my_free(vals);
      return NULL;
    }
  
    errno = 0;
    pos = strtol(tokens[pos_idx], NULL, 10);
    if(pos == 0 && errno) {
      my_warn("%s:%d: could not parse position in col %d as long int: '%s'\n", 
	    __FILE__, __LINE__, pos_idx, tokens[pos_idx]);
      my_free(tokens);
      my_free(vals);
      return NULL;
    }

    errno = 0;
    val = strtol(tokens[val_idx], NULL, 10);
    if(val == 0 && errno) {
      my_warn("%s:%d: could not parse value in col %d as long int: '%s'\n",
	      __FILE__, __LINE__, val_idx, tokens[val_idx]);
      my_free(tokens);
      my_free(vals);
      return NULL;
    }

    if(val > SHRT_MAX) {
      my_warn("%s:%d: setting value at position %ld (%ld) to maximum "
	      "allowed for int16 (%d)\n", __FILE__, __LINE__, 
	      pos, val, SHRT_MAX);
      val = SHRT_MAX;
    }

    if((pos > chr_len) || (pos < 1)) {
      my_warn("%s:%d: skipping pos %ld "
	      "(outside chromosome range 1-%ld)\n", __FILE__, __LINE__,
	      pos, chr_len);
    } else {
      vals[pos-1] = (short)val;
      count++;

      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }
    }
  }

  fprintf(stderr, "last position was %ld (%ld bp from chr end)\n", 
	  pos, chr_len - pos);
  my_free(tokens);
  fprintf(stderr, "\n");
  gzclose(gzf);
  
  return vals;
}
