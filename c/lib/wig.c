
#include <zlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

#include "wig.h"
#include "memutil.h"
#include "util.h"
#include "err.h"

/**
 * Function adapted from Michael Hoffman's genomedata wiggle parser:
 * http://noble.gs.washington.edu/proj/genomedata/
 *
 * Returns 0 on success, WIG_ERR on failure.
 *
 * Type can be either WIG_TYPE_VAR or WIG_TYPE_FIX depending on whether
 * the wiggle format is fixed or "variable"
 */
int parse_wiggle_header(char *line, char **chrom, int *type,
			long *start, long *step, long *span) {
  /* mallocs chrom; caller must free() it */
  /* start and step may be null pointers */

  char *save_ptr;
  char *token;
  char *line_no_id;
  char *loc_eq;
  char *key;
  char *val;

  line_no_id = NULL;
  *chrom = NULL;
  
  if(strncmp(WIG_ID_FIX, line, strlen(WIG_ID_FIX)) == 0) {
    *type = WIG_TYPE_FIX;
    line_no_id = &line[strlen(WIG_ID_FIX)];
  } else if (strncmp(WIG_ID_VAR, line, strlen(WIG_ID_VAR)) == 0) {
    *type = WIG_TYPE_VAR;
    line_no_id = &line[strlen(WIG_ID_VAR)];
  } else if(strncmp(WIG_ID_START, line, strlen(WIG_ID_START)) == 0) {
    /* bastard wiggle format that I used when I wrote gcons
     * LLRs to disk. Header lines just contain: start=<start>, but
     * rest is same as fixed wiggle format
     */
    *type = WIG_TYPE_FIX;
    line_no_id = line;
  } else {
    my_warn("%s:%d: unknown format on line '%s', expected '%s' or '%s'",
	    __FILE__, __LINE__, line, WIG_ID_FIX, WIG_ID_VAR);
    return WIG_ERR;
  }

  /* strip trailing newline */
  *strchr(line, '\n') = '\0';

  /* Initialize to avoid compiler warning */
  save_ptr = NULL;

  /* set to defaults */
  *span = 1;
  if(start) {
    *start = 1;
  }
  if(step) {
    *step = 1;
  }
  if(span) {
    *span = 1;
  }

  /* extra set of parentheses avoids compiler warning */
  while ((token = strtok_r(line_no_id, WIG_DELIM, &save_ptr))) {
    loc_eq = strchr(token, '=');

    /* key is allocated and must be freed in each loop iteration */
    key = util_str_ndup(token, loc_eq - token);
    assert(key);

    val = loc_eq + 1;

    if (!strcmp(key, WIG_KEY_CHROM)) {
      /* everything after the equal sign in the chromosome token */
      *chrom = util_str_dup(val);
      if(!*chrom) {
	my_warn("%s:%d: could not parse chrom from wiggle header '%s'",
		__FILE__, __LINE__, line);
	return WIG_ERR;
      }

    } else if (!strcmp(key, WIG_KEY_START)) {
      errno = 0;
      *start = strtol(val, NULL, 10);
      if(errno && *start == 0) {
	my_warn("%s:%d: start is not a valid integer '%s'",
		__FILE__, __LINE__, val);
	return WIG_ERR;
      }
    } else if (!strcmp(key, WIG_KEY_STEP)) {
      errno = 0;
      *step = strtol(val, NULL, 10);
      if(errno && *step == 0) {
	my_warn("%s:%d: step is not a valid integer '%s'",
		__FILE__, __LINE__, val);
	return WIG_ERR;
      }
    } else if (!strcmp(key, WIG_KEY_SPAN)) {
      errno = 0;
      *span = strtol(val, NULL, 10);
      if(errno && *span == 0) {
	my_warn("%s:%d: span is not a valid integer '%s'",
		__FILE__, __LINE__, val);
      }
    } else {
      my_warn("%s:%d: can't understand key: %s", __FILE__, __LINE__, key);
      my_free(key);
      return WIG_ERR;
    }

    /* free key for the next loop iteration */
    my_free(key);
    line_no_id = NULL;
  }

  return 0;
}


/**
 * Reads float32 values for an entire chromosome from a wiggle file.
 * Unspecified values are set to NAN. Returns NULL on failure.
 */
float *wig_read_float32(const char *filename, const long chr_len) {
  char line[WIG_MAX_LINE];
  gzFile gzf;
  float *vals, val;
  int type;
  char *chrom;
  long pos, step, span, count;
  
  /* allocate memory, initialize values to NaN */
  vals = my_malloc(sizeof(float) * chr_len);
  long i;
  for(i = 0; i < chr_len; i++) {
    vals[i] = NAN;
  }
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
    my_free(vals);
    return NULL;
  }

  pos = 1;
  count = 0;
  while(gzgets(gzf, line, sizeof(line))) {
    if((line[0] == 'f') || line[0] == 'v' || line[0] == 's') {
      /* parse header */
      chrom = NULL;
      if(parse_wiggle_header(line, &chrom, &type, &pos, &step, &span)) {
	my_free(vals);
	return NULL;
      }
      /* fprintf(stderr, "pos=%ld, step=%ld, span=%ld", pos, step, span); */

      if(chrom) {
	my_free(chrom);
      }
    } else {
      val = 0.0;
      if(type == WIG_TYPE_FIX) {
	/* fixed format just has values */
	val = util_parse_double(line);
      }
      else if(type == WIG_TYPE_VAR) {
	/* var format has a start position followed by a value */
	char *next;
	errno = 0;
	pos = strtol(line, &next, 10);
	if(errno && pos == 0) {
	  my_warn("%s:%d: first token in var step wig line is not "
		  "valid integer: '%s'",
		  __FILE__, __LINE__, line);
	  my_free(vals);
	  return NULL;
	}
	val = util_parse_double(next);
      }
      else {
	my_warn("%s:%d: unknown wiggle format\n", __FILE__, __LINE__);
	my_free(vals);
	return NULL;
      }

      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }

      for(i = 0; i < span; i++) {
	if((pos+i > chr_len) || (pos+i < 1)) {
	  my_warn("%s:%d: skipping pos %ld"
		  "(past chromosome end %ld)\n", __FILE__, __LINE__,
		  pos+i, chr_len);
	} else {
	  vals[pos + i - 1] = val;
	}
      }
      pos += step;
    }
  }
  fprintf(stderr, "\n");
  gzclose(gzf);
  
  return vals;
}




/**
 * Reads int16 values for an entire chromosome from a wiggle file.
 * Unspecified values are set to 0. Returns NULL on failure.
 */
short *wig_read_int16(const char *filename, const long chr_len) {
  char line[WIG_MAX_LINE];
  gzFile gzf;
  short *vals;
  char *chrom;
  long pos, step, span, count, val;
  int type;


  /* allocate memory, initialize values to 0 */
  vals = my_malloc(sizeof(short) * chr_len);
  long i;
  for(i = 0; i < chr_len; i++) {
    vals[i] = 0;
  }
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
    my_free(vals);
    return NULL;
  }

  pos = 1;
  count = 0;
  while(gzgets(gzf, line, sizeof(line))) {
    if((line[0] == 'f') || line[0] == 'v') {
      /* parse header */
      chrom = NULL;
      if(parse_wiggle_header(line, &chrom, &type, &pos, &step, &span)) {
	my_free(vals);
	return NULL;
      }
      fprintf(stderr, "pos=%ld, step=%ld, span=%ld\n", pos, step, span);

      my_free(chrom);
    } else {
      val = 0;
      if(type == WIG_TYPE_FIX) {
	/* fixed format just has values */
	val = util_parse_long(line);
      }
      else if(type == WIG_TYPE_VAR) {
	/* var format has a start position followed by a value */
	char *next;
	errno = 0;
	pos = strtol(line, &next, 10);
	if(errno && pos == 0) {
	  my_warn("%s:%d: first token in var step wig line is not "
		  "valid integer: '%s'", __FILE__, __LINE__, line);
	  my_free(vals);
	  return NULL;
	}
	val = util_parse_long(next);
      } 
      else {
	my_warn("%s:%d: unknown wiggle format\n", __FILE__, __LINE__);
	my_free(vals);
	return NULL;
      }

      if(val > SHRT_MAX) {
	my_warn("%s:%d: value %ld exceeds int16 max, setting to %d", 
		__FILE__, __LINE__, val, SHRT_MAX); 
	val = SHRT_MAX;
      }
      if(val < SHRT_MIN) {
	my_warn("%s:%d: value %ld is less than int16 min, setting to %d", 
		__FILE__, __LINE__, val, SHRT_MIN); 
	val = SHRT_MIN;
      }

      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }

      for(i = 0; i < span; i++) {
	if((pos+i > chr_len) || (pos+i < 1)) {
	  my_warn("%s:%d: skipping pos %ld"
		  "(past chromosome end %ld)\n", __FILE__, __LINE__,
		  pos+i, chr_len);
	} else {
	  vals[pos + i - 1] = val;
	}
      }
      pos += step;
    }
  }
  fprintf(stderr, "\n");
  gzclose(gzf);
  
  return vals;
}






/**
 * Reads unsigned int8 values for an entire chromosome from a wiggle file.
 * Unspecified values are set to 0. Returns NULL on failure.
 */
unsigned char *wig_read_uint8(const char *filename, const long chr_len) {
  char line[WIG_MAX_LINE];
  gzFile gzf;
  unsigned char *vals;
  char *chrom;
  long pos, step, span, count, val;
  int type;

  /* allocate memory, initialize values to 0 */
  vals = my_malloc(sizeof(unsigned char) * chr_len);
  long i;
  for(i = 0; i < chr_len; i++) {
    vals[i] = 0;
  }
  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_warn("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
    my_free(vals);
    return NULL;
  }

  pos = 1;
  count = 0;
  while(gzgets(gzf, line, sizeof(line))) {
    if((line[0] == 'f') || line[0] == 'v') {
      /* parse header */
      chrom = NULL;
      if(parse_wiggle_header(line, &chrom, &type, &pos, &step, &span)) {
	my_free(vals);
	return NULL;
      }
      fprintf(stderr, "pos=%ld, step=%ld, span=%ld\n", pos, step, span);

      my_free(chrom);
    } else {
      val = 0;
      if(type == WIG_TYPE_FIX) {
	/* fixed format just has values */
	val = util_parse_long(line);
      }
      else if(type == WIG_TYPE_VAR) {
	/* var format has a start position followed by a value */
	char *next;
	errno = 0;
	pos = strtol(line, &next, 10);
	if(errno && pos == 0) {
	  my_warn("%s:%d: first token in var step wig line is not "
		  "valid integer: '%s'", __FILE__, __LINE__, line);
	  my_free(vals);
	  return NULL;
	}
	val = util_parse_long(next);
      } 
      else {
	my_warn("%s:%d: unknown wiggle format\n", __FILE__, __LINE__);
	my_free(vals);
	return NULL;
      }

      if(val > UCHAR_MAX) {
	my_warn("%s:%d: value %ld exceeds uint8 max, setting to %d", 
		__FILE__, __LINE__, val, UCHAR_MAX); 
	val = UCHAR_MAX;
      }
      if(val < 0) {
	my_warn("%s:%d: value %ld is less than uint16 min, setting to 0", 
		__FILE__, __LINE__, val);
	val = 0;
      }

      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }

      for(i = 0; i < span; i++) {
	if((pos+i > chr_len) || (pos+i < 1)) {
	  my_warn("%s:%d: skipping pos %ld"
		  "(past chromosome end %ld)\n", __FILE__, __LINE__,
		  pos+i, chr_len);
	} else {
	  vals[pos + i - 1] = val;
	}
      }
      pos += step;
    }
  }
  fprintf(stderr, "\n");
  gzclose(gzf);
  
  return vals;
}



/**
 * Writes uint8 values in wiggle format to provided gzFile
 */
void wig_gzf_write_uint8(gzFile gzf, const unsigned char *vals,
			 const char *chr_name, const long chr_len) {
  long i;

  gzprintf(gzf, "fixedStep chrom=%s start=1 step=1\n", chr_name);
  for(i = 0; i < chr_len; i++) {
    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }
    gzprintf(gzf, "%d\n", vals[i]);
  }
  
}



/**
 * Writes gzipped uint8 values in wiggle format to the file stream
 * associated with provided file descriptor.
 */
void wig_fd_write_uint8(int fd, const unsigned char *vals,
			const char *chr_name, const long chr_len) {
  gzFile gzf;
  gzf = gzdopen(fd, "wb");
  if(gzf == NULL) {
    my_err("%s:%d:unable to open file associated with provided descriptor",
	   __FILE__, __LINE__);
  }
  wig_gzf_write_uint8(gzf, vals, chr_name, chr_len);
  fprintf(stderr, "\n");
}


/**
 * Writes int8 values for an entire chromosome to a gzipped wiggle file
 * with the provided filename.
 */
void wig_write_uint8(const char *filename, const unsigned char *vals,
		     const char *chr_name, const long chr_len) {
  gzFile gzf;
  char *out_filename;

  if(!util_has_gz_ext(filename)) {
    my_warn("%s:%d: appending '.gz' to filename", __FILE__, __LINE__);
    out_filename = util_str_concat(filename, ".gz", NULL);
  } else {
    out_filename = util_str_dup(filename);
  }
  
  gzf = util_must_gzopen(out_filename, "wb");
  fprintf(stderr, "writing to wig file '%s'\n", out_filename);
  wig_gzf_write_uint8(gzf, vals, chr_name, chr_len);
  fprintf(stderr, "\n");
  gzclose(gzf);
  my_free(out_filename);
}







/**
 * Writes float32 values for an entire chromosome to a gzipped wiggle file.
 */
void wig_write_float32(const char *filename, const float *vals,
		       const char *chr_name, const long chr_len) {
  gzFile gzf;
  long i;
  char *out_filename;

  if(!util_has_gz_ext(filename)) {
    my_warn("%s:%d: appending '.gz' to filename", __FILE__, __LINE__);
    out_filename = util_str_concat(filename, ".gz", NULL);
  } else {
    out_filename = util_str_dup(filename);
  }
  
  gzf = util_must_gzopen(out_filename, "wb");

  fprintf(stderr, "writing to wig file '%s'\n", out_filename);

  gzprintf(gzf, "fixedStep chrom=%s start=1 step=1\n", chr_name);
  for(i = 0; i < chr_len; i++) {
    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }
    gzprintf(gzf, "%.3f\n", vals[i]);
  }
  gzclose(gzf);

  fprintf(stderr, "\n");
  my_free(out_filename);
}





