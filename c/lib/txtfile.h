#ifndef __TXTFILE_H__
#define __TXTFILE_H__

#define TXTFILE_MAX_LINE 1024
#define TXTFILE_ERR -1

char *txtfile_read_int8(const char *filename, const long chr_len, 
			const int pos_idx, const int val_idx);

short *txtfile_read_int16(const char *filename, const long chr_len, 
			  const int pos_idx, const int val_idx);

#endif
