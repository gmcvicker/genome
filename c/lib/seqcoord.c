
#include <string.h>
#include <stdio.h>

#include "err.h"
#include "seqcoord.h"
#include "memutil.h"
#include "util.h"

/**
 * Copies the attributes of one sequence coordinate
 * into another
 */
void seq_coord_copy(const SeqCoord *src, SeqCoord *dst) {
  if(!memcpy(dst, src, sizeof(SeqCoord))) {
    my_err("%s:%d: sequence coordinate copy failed", __FILE__, __LINE__);
  }
  if(src->seqname) {
    dst->seqname = util_str_dup(src->seqname);
  }
  dst->score = src->score;

}


static int seq_coord_cmp_helper(SeqCoord *sc1, SeqCoord *sc2,
				 int cmp_strand, int cmp_start) {
  /* should consider comparing assembly too */
  int cmp_val;
  
  /* first need to order by sequence name */
  if(sc1->chr != NULL && sc2->chr != NULL) {
    cmp_val = strcmp(sc1->chr->name, sc2->chr->name);
    if(cmp_val != 0) {
      return cmp_val;
    }
  }

  if(sc1->seqname != NULL && sc2->seqname != NULL) {
    cmp_val = strcmp(sc1->seqname, sc2->seqname);
    if(cmp_val != 0) {
      return cmp_val;
    }
  }

  if(cmp_strand) {
    /* order with rev strand before fwd strands next*/
    if(sc1->strand < sc2->strand) {
      return -1;
    }
    if(sc1->strand > sc2->strand) {
      return 1;
    }
  }

  if(cmp_start) {
    /* sort by start position */
    if(sc1->start == sc2->start) {
      return 0;
    }
    if(sc1->start < sc2->start) {
      return -1;
    }
  } else {
    /* sort by end position */
    if(sc1->end == sc2->end) {
      return 0;
    }
    if(sc1->end < sc2->end) {
      return -1;
    }
  }

  return 1;  
}



/**
 * Comparison function for SeqCoords used for sorting. First
 * chromosome names are compared (if both names are non-NULL) and
 * coordinates are ordered alphanumerically by their seqname. Then
 * strands are compared and fwd strand coordinates are taken to be
 * "higher" than unknown strand coordinates, which are "higher" than
 * negative strand coordinates. If the seqnames and strands of the two
 * coordinates to compare are equal, the START of each of the
 * coordinates is compared.
 *
 * Currently the assembly version of the chromosomes are not taken
 * into account.
 */
int seq_coord_cmp(const void *p1, const void *p2) {
  SeqCoord *sc1, *sc2;

  sc1 = (SeqCoord *)p1;
  sc2 = (SeqCoord *)p2;

  return seq_coord_cmp_helper(sc1,sc2, TRUE, TRUE);
}



/**
 * Comparison function for SeqCoords used for sorting. First
 * chromosome names are compared (if both names are non-NULL) and
 * coordinates are ordered alphanumerically by their seqname. Then
 * strands are compared and fwd strand coordinates are taken to be
 * "higher" than unknown strand coordinates, which are "higher" than
 * negative strand coordinates. If the seqnames and strands of the two
 * coordinates to compare are equal, the END of each of the
 * coordinates is compared.
 */
int seq_coord_cmp_end(const void *p1, const void *p2) {
  SeqCoord *sc1, *sc2;

  sc1 = (SeqCoord *)p1;
  sc2 = (SeqCoord *)p2;

  return seq_coord_cmp_helper(sc1,sc2, TRUE, FALSE);
}



/**
 * Returns the cumulative length of all of the coords in
 * the provided array
 */
long seq_coord_array_len(SeqCoord *c, long num_coords) {
  long len,i;

  len = 0;
  for(i = 0; i < num_coords; i++) {
    len += seq_coord_len(&c[i]);
  }

  return len;
}


/**
 * Comparison function for SeqCoords used for sorting. First chr names
 * are compared (if both names are non-NULL) and coordinates are
 * ordered alphanumerically by their chromosome name. If the
 * chromosome names of the two coordinates are equal the start of each
 * of the coordinates is compared.
 *
 * This is the same as the seq_coord_cmp function but strand
 * information is not used.
 */
int seq_coord_cmp_nostrand(const void *p1, const void *p2) {
  SeqCoord *sc1, *sc2;

  sc1 = (SeqCoord *)p1;
  sc2 = (SeqCoord *)p2;

  return seq_coord_cmp_helper(sc1,sc2, FALSE, TRUE);  
}




void seq_coord_write(FILE *fh, SeqCoord *sc) {
  fprintf(fh, "%s:%ld-%ld(%c)", (sc->chr == NULL) ? "" : sc->chr->name,
	  sc->start, sc->end, strand_to_char(sc->strand));
}




/**
 * Returns true if the provided coordinates overlap, false otherwise.
 * If the cmp_strand argument is true, coordinates are only considered
 * overlapping if their strands match. If both coordinates have
 * non-NULL chromosomes, the chromosome name is used for comparison as
 * well.
 */
int seq_coord_ovlp(SeqCoord *sc1, SeqCoord *sc2, int cmp_strand) {
  if(sc1->chr != NULL && sc2->chr != NULL) {
    if(strcmp(sc1->chr->name, sc2->chr->name) != 0) {
      return FALSE;
    }
  } else if(sc1->seqname != NULL && sc2->seqname != NULL) {
    if(strcmp(sc1->seqname, sc2->seqname) != 0) {
      return FALSE;
    }
  }

  if(cmp_strand) {
    if(sc1->strand != sc2->strand) {
      return FALSE;
    }
  }

  if(sc1->start <= sc2->end && sc1->end >= sc2->start) {
    return TRUE;
  }

  return FALSE;
}



/**
 * Frees an array of SeqCoords. Does not free associated chromosomes.
 */
void seq_coord_array_free(SeqCoord *scs, long num) {
  int i;

  if(num > 0) {
    for(i = 0; i < num; i++) {
      if(scs[i].seqname != NULL) {
	my_free(scs[i].seqname);
      }
    }

    my_free(scs);
  }
}






/**
 * Reads sequence coordinates from a BED file. The format of the file is 
 * expected to be:
 *    chromosome start end [name [score [strand]]]
 *
 * The name attribute is not used by this function, but must be present
 * as a placeholder if the other optional attributes are to be used.
 * 
 * Other information on the BED lines, including strand, 
 * is ignored by this function.
 * 
 * BED coordinates start at 0, but we use inclusive coordinates that 
 * start at 1, so +1 is added to start values in the file.
 *
 * If a chromosome is provied then only coordinates matching that
 * chromosome are returned, otherwise all coordinates are returned. If
 * the chromosome is provided, then the returned coords have their chr
 * attribute set to point to it and their seqname set to NULL. Otherwise
 * the returned coords have chr set to NULL and seqname set to the
 * name of the chromosome they are on.
 *
 * 
 */
SeqCoord *seq_coord_read_bed(const char *filename, Chromosome *chr, 
			     long *n_coord) {
  FILE *f;
  SeqCoord *coords;
  char buf[SEQ_COORD_MAX_BED_LINE];
  char chr_name_buf[SEQ_COORD_MAX_BED_LINE];
  char name_buf[SEQ_COORD_MAX_BED_LINE];
  char strand_char;
  long i, len;
  int n_tok;

  /** TODO: should make handle gzipped files **/

  f = util_must_fopen(filename, "r");

  /* first pass, count number of coords that match chromosome name */
  *n_coord = 0;
  while(fgets(buf, sizeof(buf), f)) {
    len = strlen(buf);
    if(buf[len-1] != '\n') {
      my_err("%s:%d: BED line exceeds maximum length of %ld", __FILE__,
	     __LINE__, SEQ_COORD_MAX_BED_LINE);
    }
    if(chr) {
      if(sscanf(buf, "%s", chr_name_buf) == 0) {
	my_err("%s:%d: empty line, or first token missing", __FILE__,
	       __LINE__);
      }
      if(strcmp(chr->name, chr_name_buf) == 0) {
	*n_coord += 1;
      }
    } else {
      *n_coord += 1;
    }
  }  
  if(fseek(f, 0L, SEEK_SET) != 0) {
    my_err("%s:%d: could not rewind filehandle", __FILE__, __LINE__);
  }

  if(*n_coord == 0) {
    /* no coordinates matched chromosome */
    return NULL;
  }

  /* allocate memory for all coordinates */
  coords = my_new(SeqCoord, *n_coord);

  /* re-read file, this time parsing and setting attributes of each coordinate */
  i = 0;
  while((i < *n_coord) && (fgets(buf, sizeof(buf), f))) {
    if(i >= *n_coord) {
      break;
    }
    n_tok = sscanf(buf, "%s %ld %ld %s %lf %c", chr_name_buf,
		   &coords[i].start, &coords[i].end, name_buf,
		   &coords[i].score, &strand_char);
    
    if(n_tok < 3) {
      my_err("%s:%d: expected at least 3 tokens per BED line, got %d\n",
	     __FILE__, __LINE__, n_tok);
    }
    
    if(chr) {
      if(strcmp(chr_name_buf, chr->name)==0) {
	coords[i].chr = chr;
	coords[i].seqname = NULL;
      } else {
	/* skip--this coordinate is not on the correct chromosome */
	continue;
      }
    } else {
      coords[i].seqname = util_str_dup(chr_name_buf);
      coords[i].chr = NULL;
    }
    
    /* add one to start because BED coords start at 0 */
    coords[i].start += 1;
    
    if(n_tok < 4) {
      coords[i].score = 0.0;
    }
    if(n_tok < 5) {
      coords[i].strand = STRAND_NONE;
    } else {
      coords[i].strand = char_to_strand(strand_char);
    }
    
    i++;
  }
  if(i != *n_coord) {
    my_err("%s:%d: expected %ld matching coordinates, but only got %ld"
	   __FILE__, __LINE__, *n_coord, i);
  }  

  fclose(f);

  return coords;
}

