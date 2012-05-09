/* xbfiles - Interface with binary files. */

#ifndef __XBF_H__
#define __XBF_H__

#include <stdlib.h>
#include <stdint.h>

#define XB_MAGIC 0xCA60B175ul
#define XB_MAGIC_REV 0xB175CA60ul
#define XB_VERSION 0x0001

/* xb File structure.
 *     fixedWidthHeader
 *         magic# 		4 bytes (UNSIGNED)
 *         version              2 bytes
 *	   sizeOf(xbval_t)	2 bytes
 *         xbStranded           2 bytes (NOTE Roger said this was 2 bytes, but actually 4; strandedness, 0 no strand, 1 two strands)
 *         xbCount              2 bytes (number of xb elements + 1, We start counting from 1 here!)
 *     xbInfo                   one for each xbCount
 *         nameSize             1 byte (UNSIGNED)
 *         xbName               1*nameSize 
 *         xbSize               4 bytes (UNSIGNED) chromosome size
 *     xbRecord                 one for each xbCount (and strand if xbStranded=1)
 *         data                 sizeOf(xbval_t)*xbsize[xb_count]
 * 
 */


typedef unsigned char xbval_t;

typedef struct {
  size_t n;
  xbval_t *a;
} xbvec_t;

/* Holds header and index info from .Xb file. */
typedef struct {
  int32_t is_stranded;      /* flag: 1 if data for each strand present, 0 otherwise */
  int16_t version;	   /* version of file format */
  int16_t type_size;       /* sizeof(xbval_t); */
  int16_t count;	   /* number of chromosomes */
  char **names;            /* names of the chromosomes */
  uint32_t *sizes;         /* lengths of the chromosomes */
  xbvec_t *vec;            /* data vector */ 
  xbval_t *a;              /* The entire concatenated vector, or NULL if not used */
  void *pmmap;             /* Pointer to the mmaped area, or NULL if not mmaped */
  long long int mmap_length; /* direction to the file */
  int fd;                  /* File descriptor index */
} xblist_t;


xblist_t *xb_init(int32_t is_stranded, int16_t count, char **names, 
		  uint32_t *sizes);

void xb_free(xblist_t *xbl);

xblist_t *xb_load_mmap(const char *filename);

xblist_t *xb_init_mmap(int32_t is_stranded, int16_t count, char **names,
		       uint32_t *sizes, const char *filename);

void xb_chrom_vals(xblist_t *xbl, const char *chrom_name, xbval_t **fwd_vals,
		   xbval_t **rev_vals, long *chrom_size);

#endif
