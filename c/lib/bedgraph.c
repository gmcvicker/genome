
#include <errno.h>
#include <zlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "memutil.h"
#include "util.h"
#include "bedgraph.h"
#include "err.h"


/**
 * Reads int16 values for an entire chromosome from a tab- or
 * space-delimited text file. The column indices of the position and
 * value columns in the text file must be provided. Unspecified values
 * are set to 0. Returns NULL on error.
 */
short *bedgraph_read_int16(const char *filename, const long chr_len) {
  char buf[BEDGRAPH_MAX_LINE];
  gzFile gzf;
  short *vals;
  long start, end, i, val, count, line_num;
  char *tokens[BEDGRAPH_N_TOK], *line;

  if(chr_len < 1) {
    my_warn("%s:%d: invalid chromosmoe len (%ld)", 
	    __FILE__, __LINE__, chr_len);
    return NULL;
  }
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
    return NULL;
  }

  /* allocate memory and initialize to 0 */
  vals = (short *)my_malloc(chr_len * sizeof(short));
  memset(vals, 0, chr_len * sizeof(short));

  line_num = 0;

  /* skip past header lines looking for first line that starts with "chr" */
  line = gzgets(gzf, buf, sizeof(buf));
  while((line != NULL) && !util_str_starts_with(line, "chr")) {
    line_num++;
    line = gzgets(gzf, buf, sizeof(buf));
  }


  count = 0;
  
  while(line != NULL) {

    line_num++;
    if(util_str_split(line, tokens, BEDGRAPH_N_TOK) != BEDGRAPH_N_TOK) {
      /* error occured splitting line */
      my_warn("%s:%d: could not parse bedgraph line %d\n",
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }
  
    /* get start coordinate (starts at 0 in bedgraph)*/
    errno = 0;
    start = strtol(tokens[1], NULL, 10);
    if(start == 0 && errno) {
      my_warn("%s:%d: could not parse start from line %d\n", 
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }

    /* get end coordinate */
    errno = 0;
    end = strtol(tokens[2], NULL, 10);
    if(start == 0 && errno) {
      my_warn("%s:%d: could not parse end from line %d\n", 
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }

    /* get data value */
    errno = 0;
    val = strtol(tokens[3], NULL, 10);
    if(val == 0 && errno) {
      my_warn("%s:%d: could not parse value from line: %ld\n",
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }

    if(val > SHRT_MAX) {
      my_warn("%s:%d: setting values at positions %ld-%ld to maximum "
	      "allowed for int16 (%d)\n", __FILE__, __LINE__, 
	      start, end, val, SHRT_MAX);
      val = SHRT_MAX;
    }

    /* set values in range */
    for(i = start; i < end; i++) {
      if(i >= chr_len) {
	my_warn("%s:%d: ignoring values outside of chromosome "
		"range %ld-%ld on line %d (start=%d, end=%d)",
		__FILE__, __LINE__, 1, chr_len, line_num, start,
		end);
	break;
      }
      if(vals[i] != 0) {
	my_warn("%s:%d: value at index %d already set on line %d", __FILE__,
		__LINE__, i, line_num);
      }
      vals[i] = (short)val;


      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }
    }

    /* read next line */
    line = gzgets(gzf, buf, sizeof(buf));
  }

  fprintf(stderr, "\n");
  gzclose(gzf);
  
  return vals;
}




/**
 * Reads float32 values for an entire chromosome from a tab- or
 * space-delimited text file. The column indices of the position and
 * value columns in the text file must be provided. Unspecified values
 * are set to NAN. Returns NULL on error.
 */
float *bedgraph_read_float32(const char *filename, const long chr_len) {
  char buf[BEDGRAPH_MAX_LINE];
  gzFile gzf;
  float *vals;
  double val;
  long start, end, i, count, line_num;
  char *tokens[BEDGRAPH_N_TOK], *line;

  if(chr_len < 1) {
    my_warn("%s:%d: invalid chromosmoe len (%ld)", 
	    __FILE__, __LINE__, chr_len);
    return NULL;
  }
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, 
	    __LINE__, filename);
    return NULL;
  }

  /* allocate memory and initialize to NAN */
  vals = (float *)my_malloc(chr_len * sizeof(float));
  for(i = 0; i < chr_len; i++) {
    vals[i] = NAN;
  }


  line_num = 0;

  /* skip past header lines looking for first line that 
     starts with "chr" */
  line = gzgets(gzf, buf, sizeof(buf));
  while((line != NULL) && !util_str_starts_with(line, "chr")) {
    line_num++;
    line = gzgets(gzf, buf, sizeof(buf));
  }


  count = 0;
  
  while(line != NULL) {

    line_num++;
    if(util_str_split(line, tokens, BEDGRAPH_N_TOK) != BEDGRAPH_N_TOK) {
      /* error occured splitting line */
      my_warn("%s:%d: could not parse bedgraph line %d\n",
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }
  
    /* get start coordinate (starts at 0 in bedgraph)*/
    errno = 0;
    start = strtol(tokens[1], NULL, 10);
    if(start == 0 && errno) {
      my_warn("%s:%d: could not parse start from line %d\n", 
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }

    /* get end coordinate */
    errno = 0;
    end = strtol(tokens[2], NULL, 10);
    if(start == 0 && errno) {
      my_warn("%s:%d: could not parse end from line %d\n", 
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }

    /* get data value */
    errno = 0;
    val = strtod(tokens[3], NULL);
    if(errno) {
      my_warn("%s:%d: could not parse value from line: %ld\n",
	      __FILE__, __LINE__, line_num);
      my_free(vals);
      return NULL;
    }

    /* set values in range */
    for(i = start; i < end; i++) {
      if(i >= chr_len) {
	my_warn("%s:%d: ignoring values outside of chromosome "
		"range %ld-%ld on line %d (start=%d, end=%d)",
		__FILE__, __LINE__, 1, chr_len, line_num, start,
		end);
	break;
      }
      if(!isnan(vals[i])) {
	my_warn("%s:%d: value at index %d already set on line %d", __FILE__,
		__LINE__, i, line_num);
      }
      vals[i] = (float)val;

      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }
    }

    /* read next line */
    line = gzgets(gzf, buf, sizeof(buf));
  }

  fprintf(stderr, "done parsing bedgraph\n");
  gzclose(gzf);
  
  return vals;
}
