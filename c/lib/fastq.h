
#ifndef __FASTQ_H__

#define FASTQ_MIN_QUAL_SANGER '!'
#define FASTQ_MIN_QUAL_SOLEXA ';'
#define FASTQ_MIN_QUAL_ILLUM_1_3 '@'
#define FASTQ_MIN_QUAL_ILLUM_1_5 'B'
#define FASTQ_MIN_QUAL FASTQ_MIN_QUAL_SANGER
#define FASTQ_MAX_QUAL '~'

#define FASTQ_MAX_LINE 1024

#define FASTQ_OK 0
#define FASTQ_ERR 1
#define FASTQ_END -1

#define FASTQ_MAX_WARN 1000


typedef struct {
  char machine[FASTQ_MAX_LINE]; /* machine name */
  int run_num; /* run number */
  int lane;    /* lane number */
  int tile;    /* tile number */
  int x;       /* x-coordinate of cluster */
  int y;       /* y-coordinate of cluster */

  int type; /* read type: this is typically 1 for left read and 2 for 
	     * right read of paired reads 
	     */

  int read_len;
  char min_qual;
  char max_qual;

  char line1[FASTQ_MAX_LINE];
  char line2[FASTQ_MAX_LINE];
  char line3[FASTQ_MAX_LINE];
  char line4[FASTQ_MAX_LINE];

  int status;
} FastqSeq;



int fastq_parse_read(FastqSeq *read, gzFile f);
void fastq_report_qual_type(char min_qual, char max_qual);

#endif
