#ifndef __WIG_H__
#define __WIG_H__

#include <zlib.h>

#define WIG_MAX_LINE 1024

#define WIG_TYPE_VAR 1
#define WIG_TYPE_FIX 2
#define WIG_ID_VAR "variableStep "
#define WIG_ID_FIX "fixedStep "
#define WIG_ID_START "start="

#define WIG_DELIM " "

#define WIG_KEY_CHROM "chrom"
#define WIG_KEY_START "start"
#define WIG_KEY_STEP "step"
#define WIG_KEY_SPAN "span"

#define WIG_ERR -1

int parse_wiggle_header(char *line, char **chrom, int *type,
			long *start, long *step, long *span);

float *wig_read_float32(const char *filename, const long chr_len);

short *wig_read_int16(const char *filename, const long chr_len);

unsigned char * wig_read_uint8(const char *filename, const long chr_len);

void wig_write_uint8(const char *filename, const unsigned char *vals,
		     const char *chr_name, const long chr_len);

void wig_fd_write_uint8(int fd, const unsigned char *vals,
			const char *chr_name, const long chr_len);

void wig_gzf_write_uint8(gzFile gzf, const unsigned char *vals,
			 const char *chr_name, const long chr_len);

void wig_write_float32(const char *filename, const float *vals,
		       const char *chr_name, const long chr_len);


#endif
