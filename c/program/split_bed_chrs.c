
#include <zlib.h>
#include <ctype.h>
#include <string.h>

#include <glib.h>

#include "memutil.h"
#include "util.h"
#include "err.h"

#define BED_MAX_LINE 1000000
#define MAX_TOK MAX_LINE
#define CHR_MAX 128

char *get_output_path(char *old_path, char *chr) {
  int i;  
  char *dir, *filename, *new_filename;

  dir = util_str_dup(old_path);
    
  /* search for final '/' in path */
  i = strlen(old_path)-1;
  while(i > 0) {
    if(old_path[i] == '/') {
      break;
    }
    i--;
  }
  filename = &old_path[i+1];
  
  if(i > 0) {
    dir[i] = '\0';
  } else {
    dir[0] = '\0';
  }

  new_filename = util_str_concat(dir, "/", chr, ".bed", NULL);

  my_free(dir);

  return new_filename;
}


void close_file(void *key, void *val, void *not_used) {
  fprintf(stderr, "closing file for %s\n", (char *)key);
  fclose(val);
}

void split_bed_chrs(char **filenames, int n_filenames) {
  char line[BED_MAX_LINE], *out_filename;
  gzFile *gzf;
  FILE *out_f;
  char chr[CHR_MAX], cur_chr[CHR_MAX];
  int i, j, is_new_chr;
  long count;
  GHashTable *file_tab;

  
  file_tab = g_hash_table_new_full(g_str_hash, g_str_equal,
				   NULL, NULL);

  out_filename = NULL;
  out_f = NULL;

  for(j = 0; j < CHR_MAX; j++) {
    cur_chr[j] = '\0';
  }

  for(i = 0; i < n_filenames; i++) {
    fprintf(stderr, "\nprocessing file %s\n", filenames[i]);

    gzf = gzopen(filenames[i], "rb");
    if(!gzf) {
      my_err("%s:%d: could not open file %s", __FILE__, __LINE__, 
	     filenames[i]);
    }

    count = 0;

    while(gzgets(gzf, line, sizeof(line))) {
      /* read chromosome name */
      for(j = 0; j < CHR_MAX; j++) {
	if(line[j] == '\0' || isspace(line[j])) {
	  break;
	}
	chr[j] = line[j];
      }
      chr[j] = '\0';

      is_new_chr = FALSE;
      if(strcmp(chr, cur_chr) != 0) {
	/* this chr differs from last one we've seen */
	is_new_chr = TRUE;
	strncpy(cur_chr, chr, CHR_MAX);
      }

      if(is_new_chr) {
	if(strncmp(chr, "chr", 3) != 0) {
	  fprintf(stderr, "WARNING: strange chromosome "
		  "name LINE:\n%s\n", line);
	}

	/* get output file for chromosome */
	if(g_hash_table_lookup(file_tab, chr)) {
	  /* we've already started writing to this file */
	  fprintf(stderr, "%s\n", chr);
	  out_f = g_hash_table_lookup(file_tab, chr);
	} else {
	  /* open a new file */
	  out_filename = get_output_path(filenames[0], chr);
	  out_f = fopen(out_filename, "w");
	  fprintf(stderr, "%s\n", chr);
	  if(!out_f) {
	    my_err("%s:%d: could not open output file %s", __FILE__, __LINE__,
		   out_filename);
	  }
	  g_hash_table_insert(file_tab, util_str_dup(chr), out_f);
	}
      }

      /* write line to file */
      /* gzprintf(out_f, "%s", line); */
      fprintf(out_f, "%s", line);
	 
      /* fprintf(stderr, "%s", line); */

      count++;
      if(count > 10000) {
	fprintf(stderr, ".");
	count = 0;
      }
    }
  }
  fprintf(stderr, "\n");

  g_hash_table_foreach(file_tab, close_file, NULL);
  
  return;
}



int main(int argc, char **argv) {
  if(argc < 2) {
    fprintf(stderr, "usage: %s <file1.bed> [<file2.bed> ...]\n", argv[0]);
    exit(2);
  }

  split_bed_chrs(&argv[1], argc-1);
  
  return 0;
}


