
#include <stdint.h>
#include <sys/mman.h>
#include <string.h>

#include "err.h"
#include "xbf.h"
#include "memutil.h"
#include "util.h"

/**
 * initializes memory for xb data consisting of count chromosomes with
 * specified names and sizes
 */
xblist_t *xb_init(int32_t is_stranded, int16_t count, char **names, 
		  uint32_t *sizes) {
  xblist_t *xbl;
  int j, k;

  xbl = my_new(xblist_t, 1);
  xbl->version = XB_VERSION;
  xbl->type_size = sizeof(xbval_t);
  xbl->count = count;
  xbl->vec = my_new(xbvec_t,  xbl->count << is_stranded);
  xbl->names = my_new(char *, xbl->count);
  xbl->sizes = my_new(uint32_t, xbl->count);

  fprintf(stderr, "before loop %d\n", count);
  for(j = 0; j < xbl->count; j++) {
    fprintf(stderr, "%d, %s %d\t", j, names[j], sizes[j]);
    xbl->names[j] = util_str_dup(names[j]);
    xbl->sizes[j] = sizes[j];
    fprintf(stderr, xbl->names[j], xbl->sizes[j]);
  }

  if(is_stranded){
    for(j = 0; j < xbl->count; j++){ 
      k = j << 1;
      /* allocate space for forward-strand vector of values */
      xbl->vec[k].n = sizes[k]; 
      xbl->vec[k].a = my_new(xbval_t, sizes[k]);

      /* allocate space for reverse-strand vector of values */
      xbl->vec[k+1].n = sizes[k+1]; 
      xbl->vec[k+1].a = my_new(xbval_t, sizes[k+1]);
    }
  } else {
    for(j = 0; j < xbl->count; j++) {
      /* allocate space for unstranded vector of values */
      xbl->vec[j].n = sizes[j]; 
      xbl->vec[j].a = my_new(xbval_t, sizes[j]);
    }
  }
  xbl->a = NULL;
  xbl->pmmap = NULL;
  xbl->mmap_length = 0;
  return xbl;
}


/**
 * Frees memory in xblist datastructure and unmaps memory-mapped file
 */
void xb_free(xblist_t *xbl) {
  int j;
  
  if(xbl->pmmap == NULL){ 
    for(j=0; j < (xbl->count << xbl->is_stranded); j++) {
      my_free(xbl->vec[j].a);
    }
  }
  else{
    /* free memory-mapped file */
    munmap(xbl->pmmap, xbl->mmap_length);
  }
  my_free(xbl->vec);
  my_free(xbl->sizes);
  my_free(xbl->names);
  
  my_free(xbl);
}



/**
 * Given the xblist sets pointers to values for the specifed chromosome
 * specified chromosome region. If the data is single-stranded the
 * reverse ptr is set to NULL.
 */
void xb_chrom_vals(xblist_t *xbl, const char *chrom_name, xbval_t **fwd_vals,
		   xbval_t **rev_vals, long *chrom_size) {
  int i;

  for(i = 0; i < xbl->count; i++) {
    if(strcmp(xbl->names[i], chrom_name) == 0) {
      break;
    }
  }

  if(i == xbl->count) {
    /* did not find a matching chromsomoe */   
    *fwd_vals = NULL;
    *rev_vals = NULL;
    *chrom_size = 0;
    my_err("%s:%d: chromosome '%s' is not present", __FILE__, __LINE__,
	   chrom_name);
  }

  /* set chromosome size */
  *chrom_size = xbl->sizes[i];

  if(xbl->is_stranded) {
    /* forward and reverse values are offset */
    i = i * 2;
    /* set ptr to forward and reverse values */
    *fwd_vals = xbl->vec[i].a;
    *rev_vals = xbl->vec[i+1].a;
  } else {
    *fwd_vals = xbl->vec[i].a;
    *rev_vals = NULL;
  }

  return;
}



/** 
 * Reads header and memory-maps contents of binary xb file. Returns the
 * xblist file.
 */
xblist_t *xb_load_mmap(const char *filename) {
  FILE *f;
  xblist_t *xbl;
  uint32_t magic;
  uint8_t name_size;
  unsigned long long int total_size = 0;
  long long int offset = 0;
  int j;
  int fdi;

  xbl = my_new(xblist_t, 1);
  f = util_must_fopen(filename, "rb");

  util_fread_one(f, magic);

  fprintf(stderr, "magic: %u (expected: %lu)\n", magic, XB_MAGIC);

  /* verify magic number is correct */
  if(magic != XB_MAGIC) {
    if(magic == XB_MAGIC_REV) {
      my_err("%s:%d: system appears to be big endian and byte-swapping "
	     "not yet implemented", __FILE__, __LINE__);
    } else {
      my_err("%s:%d: magic number mismatch, wrong file type?",
	     __FILE__, __LINE__);
    }
  }

  util_fread_one(f, xbl->version);
  util_fread_one(f, xbl->type_size);
  util_fread_one(f, xbl->is_stranded);
  util_fread_one(f, xbl->count);

  fprintf(stderr, "version: %d, type_size: %d, is_stranded: %d, count: %u\n",
	  xbl->version, xbl->type_size, xbl->is_stranded, xbl->count);
  
  fprintf(stderr, "there are %d chromosomes\n", xbl->count);
  xbl->vec = my_new(xbvec_t, xbl->count << xbl->is_stranded);
  xbl->names = my_new(char *, xbl->count);
  xbl->sizes = my_new(uint32_t, xbl->count);
    
  /* sequence specific info and calculating offsets...*/
  for(j = 0; j < xbl->count; j++){ //chr 0 reserved...

    /* read size of chromosome name */
    util_fread_one(f, name_size);

    /* read chromosome name */
    xbl->names[j] = my_new(char, name_size+1);
    util_must_fread(f, xbl->names[j], name_size);
    xbl->names[j][name_size] = '\0';

    /* read length of chromosome */
    util_fread_one(f, xbl->sizes[j]);

    total_size += xbl->sizes[j];
  }

  offset = ftell(f);

  /* calculate total size of file */
  xbl->mmap_length = offset + 
    (((unsigned long long)total_size << (xbl->is_stranded))
     * (unsigned long long)xbl->type_size);

  /* now memory-map the file */
  fprintf(stderr, "memory mapping %lld bytes\n", xbl->mmap_length);
  rewind(f);
  fdi = fileno(f);

  /* if read/write permission is needed, then need to change to this:
   * xbl->pmmap = mmap(0, xbl->mmap_length, PROT_READ|PROT_WRITE,
   * MAP_SHARED, fdi, 0);
   */

  xbl->pmmap = mmap(0, xbl->mmap_length, PROT_READ, MAP_SHARED, fdi, 0);

  if(xbl->pmmap == MAP_FAILED) {
    perror("mmap");
    my_err("%s:%d: memory mapping of file failed\n", __FILE__, __LINE__);
  }

  /* a is a pointer to the start of all of the data vectors */
  xbl->a = (xbl->pmmap + offset);

  if(xbl->is_stranded) {
    total_size = 0;
    for(j = 0; j < xbl->count; j++){ //chr 0 reserved...    
      int k = j << 1;
      xbl->vec[k].a = &(xbl->a[total_size]);
      xbl->vec[k].n = xbl->sizes[j];
      total_size += (unsigned long long)xbl->sizes[j];
      xbl->vec[k+1].a = &(xbl->a[total_size]);
      xbl->vec[k+1].n = xbl->sizes[j];
      total_size += (unsigned long long) xbl->sizes[j];
    }
  } else {
    total_size=0;
    for (j = 0; j < xbl->count; j++){ //chr 0 reserved...    
      xbl->vec[j].a = &(xbl->a[total_size]);
      xbl->vec[j].n = xbl->sizes[j];
      total_size += (unsigned long long) xbl->sizes[j];
    }
  }

  /* fclose(f); */
  
  return xbl;
}



/**
 * Creates a new memory-mapped file and writes its header. Returns a pointer to
 * the memory-mapped xblist structure
 */
xblist_t *xb_init_mmap(int32_t is_stranded, int16_t count, char **names,
		       uint32_t *sizes, const char *filename) {
  xbval_t val;
  FILE *f;
  uint32_t magic = XB_MAGIC;
  int16_t version = XB_VERSION;
  int16_t type_size;
  uint8_t name_size;
  unsigned long long total_size = 0;
  int j;  
  
  type_size = sizeof(xbval_t);

  f = fopen(filename,"wb");
  fprintf(stderr, "initializing header of %s...\n", filename);

  /* Write out fixed parts of header. */
  util_fwrite_one(f, magic);
  util_fwrite_one(f, version);
  util_fwrite_one(f, type_size);
  util_fwrite_one(f, is_stranded);
  util_fwrite_one(f, count);
  
  /* Sequence specific info and calculating offsets...*/
  for(j = 0; j < count; j++) { //chr 0 reserved...
    name_size = strlen(names[j]);
    util_must_fwrite(f, names[j], name_size);
    util_fwrite_one(f, sizes[j]);
    total_size += (unsigned long long) sizes[j];
  }
  fprintf(stderr, "..completed\n");

  fprintf(stderr, "initializing data space for %lld bytes...\n",
	  (total_size << is_stranded) * type_size);

  /* seek to what would be end of file (-1 element), and write single element */
  val = 0;
  fseek(f, ((total_size << is_stranded) - 1) * type_size, SEEK_CUR);
  fwrite(&val,type_size, 1, f);
  fclose(f);

  fprintf(stderr, "..completed\n");
  
  return xb_load_mmap(filename);
}


