
#include <zlib.h>
#include <ctype.h>
#include <string.h>

#include "memutil.h"
#include "util.h"
#include "err.h"
#include "wig.h"

#define MAX_LINE 1024
#define MAX_TOK MAX_LINE


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

  if(util_has_gz_ext(filename)) {
    if(strlen(dir) > 0) {
      new_filename = util_str_concat(dir, "/", chr, "_", filename, NULL);
    } else {
      new_filename = util_str_concat(chr, "_", filename, NULL);
    }
  } else {
    /* add .gz extension to filename */
    if(strlen(dir) > 0) {
      new_filename = util_str_concat(dir, "/", chr, "_", 
				     filename, ".gz", NULL);
    } else {
      new_filename = util_str_concat(chr, "_", filename, ".gz", NULL);
    }
  }

  my_free(dir);

  return new_filename;
}


void split_wig_chrs(char *filename) {
  char line[WIG_MAX_LINE], *out_filename, *header;
  gzFile gzf, out_gzf;
  int type;
  char *chr, *cur_chr;
  long pos, step, span, count;

  
  gzf = gzopen(filename, "rb");
  if(!gzf) {
    my_err("%s:%d: could not open file %s", __FILE__, __LINE__, filename);
  }

  out_filename = NULL;
  cur_chr = NULL;
  out_gzf = NULL;

  pos = 1;
  count = 0;
  while(gzgets(gzf, line, sizeof(line))) {
    if(line[0] == 'f' || line[0] == 'v') {
      /* parse header */
      header = util_str_dup(line);
      chr = NULL;
      if(parse_wiggle_header(line, &chr, &type, &pos, &step, &span)) {
	my_err("%s:%d: failed to parse wiggle header. line:\n%s",
	       __FILE__, __LINE__, line);
      }

      if((cur_chr == NULL) || (strcmp(cur_chr, chr) != 0)) {
	/* starting a new chromosome */
	if(cur_chr) {
	  my_free(cur_chr);
	  /* close old output file */
	  gzclose(out_gzf);
	  my_free(out_filename);
	}
	cur_chr = chr;
	fprintf(stderr, "%s\n", cur_chr);

	/* open new file */
	out_filename = get_output_path(filename, cur_chr);
	out_gzf = gzopen(out_filename, "wb");
	if(!out_gzf) {
	  my_err("%s:%d: could not open output file %s", __FILE__, __LINE__,
		out_filename);
	}
      }
      /* write header line to file */
      gzprintf(out_gzf, "%s", header);
      my_free(header);
    } else {
      gzprintf(out_gzf, line);

      count++;
      if(count > 1000000) {
	fprintf(stderr, ".");
	count = 0;
      }
    }
  }
  fprintf(stderr, "\n");
  gzclose(gzf);

  if(cur_chr) {
    my_free(cur_chr);
    /* close old output file */
    gzclose(out_gzf);
    my_free(out_filename);
  }
  
  return;
}





int main(int argc, char **argv) {
  char *in_filename;

  if(argc != 2) {
    fprintf(stderr, "usage: %s <filename>\n", argv[0]);
    exit(2);
  }

  in_filename = argv[1];
  split_wig_chrs(in_filename);
  
  return 0;
}


