#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "fastq.h"
#include "err.h"
#include "util.h"

static int warn_count = 0;



void fastq_check_line_len(FastqSeq *read, char *line, gzFile f) {
  size_t len, n;
  char c;

  len = strlen(line);
  if(len == 0) {
    return;
  }
  if(line[len-1] == '\n') {
    return;
  }
  
  /* line did not terminate with a '\n' */
  my_warn("%s:%d: line did not terminate with '\\n':  \n'%s'\n", 
          __FILE__, __LINE__, line, len);

  read->status = FASTQ_ERR;

  /* seek in file until next '\n' is found */
  n = 0;
  while((c = gzgetc(f)) != -1) {
    n++;

    if(n < 10) {
      if(isprint(c)) {
        fprintf(stderr, "  extra character %ld: '%c'\n", n, c);
      }
      else {
        fprintf(stderr, "  unprintable extra character %ld: '\\%d'\n", n, c);
      }
    } else if(n == 10) {
      fprintf(stderr, "  ...\n");      
    }
    
    if(c == '\n') {
      fprintf(stderr, "  read %ld extra characters to reach end of line\n", n);
      return;
    }
  }

  fprintf(stderr, "  read %ld extra characters to reach end of file\n", n);
  return;
}


static void seek_next_header(gzFile f) {
  char c1, c2;
  size_t n;

  c1 = '\n';
  c2 = gzgetc(f);

  /* search for the next line that starts with '@' */
  n = 0;
  while(c2 != -1) {
    if((c1 == '\n') && (c2 == '@')) {
      /* backup one byte to put '@' back on file stream */
      /* gzseek(f, -1, SEEK_CUR);*/
      gzungetc(c2, f);

      if(warn_count < FASTQ_MAX_WARN) {
	fprintf(stderr, "skipped %ld bytes to find next fastq header line\n", 
		n);
      }
      return;
    }
    c1 = c2;
    c2 = gzgetc(f);
    n++;
  }
  fprintf(stderr, "skipped %ld bytes at end of file\n", n);
}


/** 
 * Reads the four lines of the fastq record
 */
int fastq_read_lines(FastqSeq *read, gzFile f) {
  /* read the four lines that make up fastq record */
  if(gzgets(f, read->line1, FASTQ_MAX_LINE) == NULL) {
    /* end of file */
    read->status = FASTQ_END;
    read->line1[0] = '\0';
    read->line2[0] = '\0';
    read->line3[0] = '\0';
    read->line4[0] = '\0';
    return FASTQ_END;
  }

  /* check that this line was a header starting with '@' */
  if(read->line1[0] != '@') {
    if(warn_count < FASTQ_MAX_WARN) {
      warn_count += 1;
      my_warn("%s:%d: fastq header line does not start with '@'",
	      __FILE__, __LINE__);
    }
    read->status = FASTQ_ERR;
    read->line2[0] = '\0';
    read->line3[0] = '\0';
    read->line4[0] = '\0';
    /* move ahead in file to next line that starts with '@' */
    seek_next_header(f);
    return read->status;
  }

  fastq_check_line_len(read, read->line1, f);
  util_str_rstrip(read->line1);

  /* read second line */
  if(gzgets(f, read->line2, FASTQ_MAX_LINE) == NULL) {
    /* end of file */
    my_warn("%s:%d: fastq file ended mid-record\n",
	    __FILE__, __LINE__);
    read->status = FASTQ_ERR;
    read->line2[0] = '\0';
    read->line3[0] = '\0';
    read->line4[0] = '\0';
    return FASTQ_ERR;
  }
  fastq_check_line_len(read, read->line2, f);
  util_str_rstrip(read->line2);

  /* read third line */
  if(gzgets(f, read->line3, FASTQ_MAX_LINE) == NULL) {
    /* end of file */
    my_warn("%s:%d: fastq file ended mid-record\n",
	    __FILE__, __LINE__);
    read->status = FASTQ_ERR;
    read->line3[0] = '\0';
    read->line4[0] = '\0';
    return FASTQ_ERR;
  }
  fastq_check_line_len(read, read->line3, f);
  util_str_rstrip(read->line3);

  /* read fourth line */
  if(gzgets(f, read->line4, FASTQ_MAX_LINE) == NULL) {
    /* end of file */
    my_warn("%s:%d: fastq file ended mid-record\n", __FILE__, __LINE__);
    read->status = FASTQ_ERR;
    read->line4[0] = '\0';
    return FASTQ_ERR;
  }
  fastq_check_line_len(read, read->line4, f);
  util_str_rstrip(read->line4);

  return read->status;
}


/**
 * Checks that the header has expected 7 fields. This could be changed
 * to allow for variety of header types.
 */
static int check_header(FastqSeq *read) {
  char *offset;
  size_t assigned;

  /* parse read attributes from header line, assuming it has standard
   * formatting. Example header line:
   * IPAR1:1:2:18330:12837#0/1
   */
  offset = index(read->line1, ':');
  if(offset == NULL) {
    assigned = 0;
  } else {
    /* replace first ':' with ' ', so that string directive of
     * sscanf stops after parsing machine name
     */
    offset[0] = ' ';
    /* parse attributes */
    assigned = sscanf(read->line1, "@%s %d:%d:%d:%d#%d/%d",
		      read->machine, &read->lane, &read->tile, &read->x, 
		      &read->y,  &read->run_num, &read->type);
    offset[0] = ':';
  }
  if(assigned != 7) {
    /* failed to completely parse header */
    my_warn("%s:%d: could only parse %d out of 7 expected fields from header",
	    __FILE__, __LINE__, assigned);
    read->status = FASTQ_ERR;
  }

  return read->status;
}


int fastq_check_seq(FastqSeq *read) {
  int i, j;
  char c;
  int err;
  /* string of valid nucleotide identifiers, including ambiguity codes */
  static const char *valid_nucs = "ATCGNatcgnMRWSYKmrwsyk";

  err = FALSE;
  for(i = 0; i < read->read_len; i++) {
    c = read->line2[i];

    j = 0;
    while(c != valid_nucs[j]) {
      if(valid_nucs[j] == '\0') {
	if(warn_count < FASTQ_MAX_WARN) {
	  warn_count += 1;
	  my_warn("%s:%d: read contains invalid base '%c'", 
		  __FILE__, __LINE__, c);
	}
	err = TRUE;
	break;
      }
      j++;
    }

    if(err) {
      read->status = FASTQ_ERR;
      break;
    }
  }

  return read->status;
}



/**
 * Checks that quality characters fall within valid range
 */
int fastq_check_qual(FastqSeq *read) {
  int i;
  char c;

  read->min_qual = -1;
  read->max_qual = -1;

  for(i = 0; i < read->read_len; i++) {
    c = read->line4[i];

    if(read->min_qual == -1) {
      read->min_qual = c;
      read->max_qual = c;
    } else {
      if(c < read->min_qual) {
	read->min_qual = c;
      }
      if(c > read->max_qual) {
	read->max_qual = c;
      }
    }
  }

  if(read->min_qual < FASTQ_MIN_QUAL) {
    my_warn("%s:%d: read has invalid quality value with ascii code %d",
	    __FILE__, __LINE__, read->min_qual);
    read->status = FASTQ_ERR;
  }
  if(read->max_qual > FASTQ_MAX_QUAL) {
    my_warn("%s:%d: read has invalid quality value with ascii code %d",
	    __FILE__, __LINE__, read->max_qual);
    read->status = FASTQ_ERR;
  }

  return read->status;
}


/**
 * Parses a read in fastq format.
 * Returns FASTQ_END at end of file, FASTQ_OK on success, FASTQ_ERR on problem
 */
int fastq_parse_read(FastqSeq *read, gzFile f) {
  size_t qual_len;

  read->status = FASTQ_OK;

  fastq_read_lines(read, f);
  if(read->status != FASTQ_OK) {
    return read->status;
  }
  
  /* check_header(read); */
  /* if(read->status != FASTQ_OK) { */
  /*   return read->status; */
  /* } */

  /* third line should start with '+' separator */
  if(read->line3[0] != '+') {
    if(warn_count < FASTQ_MAX_WARN) {
      warn_count += 1;
      my_warn("%s:%d: third line does not start with '+'", 
	      __FILE__, __LINE__);
    }
    read->status = FASTQ_ERR;
    return read->status;
  }
   
  /* check length of read and quality */
  read->read_len = strlen(read->line2);
  qual_len = strlen(read->line4);
  
  if(read->read_len < 1) {
    if(warn_count < FASTQ_MAX_WARN) {
      warn_count += 1;
      my_warn("%s:%d: read has no bases\n", __FILE__, __LINE__);
    }
    return read->status;
  }

  /* next line should be quality scores */
  if(read->read_len != qual_len) {
    if(warn_count < FASTQ_MAX_WARN) {
      warn_count += 1;
      my_warn("%s:%d: read len (%ld) does not match quality score len (%ld)",
	      __FILE__, __LINE__, read->read_len, qual_len);
    }
    read->status = FASTQ_ERR;
    return read->status;
  }

  fastq_check_seq(read);
  fastq_check_qual(read);

  return read->status;
}



void fastq_report_qual_type(char min_qual, char max_qual) {
  fprintf(stderr, "\n");
  fprintf(stderr, "guessing quality format:\n");

  if((min_qual == -1) || (max_qual == -1)) {
    fprintf(stderr, "  no valid quality scores to guess quality type from\n");
  }

  fprintf(stderr, "  min_qual:%c, max_qual:%c\n", min_qual, max_qual);
  if(min_qual < FASTQ_MIN_QUAL_SOLEXA) {
    fprintf(stderr, "  quality vals appear to be Sanger / Illum 1.8+ format"
	    " (Phred+33)\n");
    
    if(max_qual >= 'h') {
      my_warn("%s:%d: quality vals may be mix of Phred+33 and Phred+64\n"
	      "         You should probably fix this.", __FILE__, __LINE__);
    }
  } else {
    if(min_qual < FASTQ_MIN_QUAL_ILLUM_1_3) {
      my_warn("%s:%d: quality vals appear to be OLD solexa format, "
	      "may need to convert prior to processing.\n", __FILE__, __LINE__);
    }
    else if(min_qual < FASTQ_MIN_QUAL_ILLUM_1_5) {
      fprintf(stderr, "  quality vals appear to be Illumina 1.3+ format "
	      "(Phred+64)\n  should probably use -I flag for bwa aln "
	      "(relevant only if using -q argument)\n");
    }
    else {
      fprintf(stderr, "  quality vals appear to be Illumina 1.5+ format "
	      "(Phred+64)\n  should probably use -I flag for bwa aln "
	      "(relevant only if using -q argument)\n");
      
    }
  }
}



