#ifndef __BEDGRAPH_H__
#define __BEDGRAPH_H__

#define BEDGRAPH_MAX_LINE 1024
#define BEDGRAPH_N_TOK 4


short *bedgraph_read_int16(const char *filename, const long chr_len);
float *bedgraph_read_float32(const char *filename, const long chr_len);


#endif
