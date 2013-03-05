

#include <glib.h>
#include <limits.h>
#include <string.h>

#include <zlib.h>

#include "seq.h"
#include "memutil.h"
#include "nuc.h"
#include "util.h"

#define SEQ_DIR "/mnt/lustre/home/gmcvicker/data/seq/hg18_nb"
#define SNP_DIR "/data/share/10_IND/IMPUTE/hg18"

#define UNMAPPABLE 255

void add_kmer_count(GHashTable *count_tab, char *kmer_str,
		    unsigned char *nucs, unsigned char *ref_nucs,
		    unsigned char *alt_nucs, int kmer_size, int offset,
		    int depth) {
  int i;
  unsigned char *count_ptr;

  for(i = offset; i < kmer_size; i++) {
    kmer_str[i] = nuc_id_to_char(nucs[i]);
    
    if((ref_nucs[i] != NUC_N) && (ref_nucs[i] != NUC_GAP)) {
      /* Count versions of kmer with reference and alternate
       * nucleotide at this position. Do this recursively so
       * that all possible combinations can be considered.
       */
      kmer_str[i] = nuc_id_to_char(ref_nucs[i]);
      add_kmer_count(count_tab, kmer_str, nucs, ref_nucs, alt_nucs,
		     kmer_size, i+1, depth+1);

      /* alternate version */
      kmer_str[i] = nuc_id_to_char(alt_nucs[i]);
      add_kmer_count(count_tab, kmer_str, nucs, ref_nucs, alt_nucs,
		     kmer_size, i+1, depth+1);

      return;
    }
  }

  count_ptr = g_hash_table_lookup(count_tab, kmer_str);
  if(count_ptr == NULL) {
    /* first time we have seen this kmer */
    count_ptr = my_new(unsigned char, 1);
    *count_ptr = 1;
    g_hash_table_insert(count_tab, util_str_dup(kmer_str),
			count_ptr);
  } else {
    /* increment count for this kmer */
    if(*count_ptr < UCHAR_MAX) {
      *count_ptr += 1;
    }
  }
}


char *get_chrom_name(char *path) {
  int len, i, last_dot;
  char *chrom_name;
  
  len = strlen(path);
  chrom_name = g_new(char, len+1);

  last_dot = len-1;
  for(i=len-1; i >= 0; i--) {
    if(path[i] == '.') {
      last_dot = i;
    }
    if(path[i] == '/') {
      i = i+1;
      break;
    }
  }

  if(i < 0) {
    i = 0;
  }

  strncpy(chrom_name, &path[i], last_dot - i);
  
  return chrom_name;
}



int lookup_kmer_count(GHashTable *count_tab, char *kmer_str,
		      unsigned char *nucs, unsigned char *ref_nucs,
		      unsigned char *alt_nucs, int kmer_size, int offset,
		      int depth) {
  int i;
  unsigned char *count_ptr, ref_count, alt_count;

  for(i = offset; i < kmer_size; i++) {
    kmer_str[i] = nuc_id_to_char(nucs[i]);
    
    if(nucs[i] == NUC_N) {
      /* flag genomic regions with Ns as unmappable */
      return UNMAPPABLE;
    }

    if(ref_nucs[i] != NUC_N) {
      /* there is a SNP at this site */

      if(ref_nucs[i] == NUC_GAP) {
	/* variant is an INDEL, flag as unmappable */
	return UNMAPPABLE;
      }
      
      if((ref_nucs[i] != nucs[i]) && (alt_nucs[i] != nucs[i])) {
	/* neither reference nor alt match reference */
	/* flag as unmappable since this looks suspicious */
	return UNMAPPABLE;
      }
            
      /* get counts for versions of kmer with reference and alternate
       * nucleotide at this position. Do this recursively so
       * that all possible combinations can be considered.
       */
      kmer_str[i] = nuc_id_to_char(ref_nucs[i]);
      ref_count = lookup_kmer_count(count_tab, kmer_str, nucs, 
				    ref_nucs, alt_nucs,
				    kmer_size, i+1, depth+1);

      /* alternate version */
      kmer_str[i] = nuc_id_to_char(alt_nucs[i]);
      alt_count = lookup_kmer_count(count_tab, kmer_str, nucs,
				    ref_nucs, alt_nucs,
				    kmer_size, i+1, depth+1);

      /* return whichever count is higher */

      return (ref_count > alt_count) ? ref_count : alt_count;
    }
  }

  /* lookup the count for this kmer and return it */
  count_ptr = g_hash_table_lookup(count_tab, kmer_str);
  if(count_ptr == NULL) {
    my_err("kmer not found in lookup table: %s\n", kmer_str);
  }

  return *count_ptr;
}




void count_kmers(GHashTable *count_tab, Seq *seq,
		 unsigned char *ref_alleles, unsigned char *alt_alleles,
		 int kmer_size) {
  long i;
  char *kmer_str;
  
  /* allocate memory to hold kmer */
  kmer_str = my_new(char, kmer_size + 1);
  kmer_str[kmer_size] = '\0';

  for(i = 0; i < (seq->len - kmer_size + 1); i++) {
    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }

    add_kmer_count(count_tab, kmer_str, &seq->sym[i],
		   &ref_alleles[i], &alt_alleles[i], kmer_size, 0, 0);
  }

  my_free(kmer_str);
    
  fprintf(stderr, "\n");
}



void report_kmer_hits(gzFile gzf, GHashTable *count_tab, Seq *seq,
		      unsigned char *ref_alleles, unsigned char *alt_alleles,
		      int kmer_size) {
  long i;
  char *kmer_str;
  unsigned char count;
  
  /* allocate memory to hold kmer */
  kmer_str = my_new(char, kmer_size + 1);
  kmer_str[kmer_size] = '\0';

  for(i = 0; i < (seq->len - kmer_size + 1); i++) {
    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }

    count = lookup_kmer_count(count_tab, kmer_str, &seq->sym[i],
			      &ref_alleles[i], &alt_alleles[i], kmer_size, 
			      0, 0);

    gzprintf(gzf,  "%d\n", count);
  }
  
  my_free(kmer_str);
}



void read_seq(Seq *seq, char *filename) {
  gzFile gzf;
  
  fprintf(stderr, "reading sequence from %s\n", filename);
  
  gzf = util_must_gzopen(filename, "rb");
  seq_read_fasta_record(seq, gzf);
  gzclose(gzf);
}




void read_snps(char *chrom_name, unsigned char *ref_alleles, 
	       unsigned char *alt_alleles, long chrom_len) {
  char snp_filename[1024];
  char snp_line[1000000];
  char allele1[10], allele2[10];
  unsigned char ref_nuc_id, alt_nuc_id;
  gzFile gzf;
  int n_tok;
  long pos, i;

  /* initialize allele arrays */
  for(i = 0; i < chrom_len; i++) {
    ref_alleles[i] = NUC_N;
    alt_alleles[i] = NUC_N;
  }
  
  snprintf(snp_filename, sizeof(snp_filename), 
	   "%s/%s.hg18.impute2.gz", SNP_DIR, chrom_name);

  if(!util_file_exists(snp_filename)) {
    fprintf(stderr, "no SNP file for %s\n", chrom_name);
    return;
  }
  fprintf(stderr, "reading SNPs from %s\n", snp_filename);

  gzf = util_must_gzopen(snp_filename, "rb");
    
  while(gzgets(gzf, snp_line, sizeof(snp_line)) != NULL) {
    n_tok = sscanf(snp_line, "%*s %*s %ld %10s %10s",  &pos, 
		   allele1, allele2);
    
    if(n_tok != 3) {
      my_err("expected 5 tokens but got %d", n_tok);
    }
    
    if((strlen(allele1) != 1) || strlen(allele2) != 1) {
      /* this is an indel, flag it as such */
      ref_nuc_id = NUC_GAP;
      alt_nuc_id = NUC_GAP;
    } else {
      ref_nuc_id = nuc_char_to_id(allele1[0]);
      alt_nuc_id = nuc_char_to_id(allele2[0]);
    }

    if((pos < 1) || (pos > chrom_len)) {
      my_err("snp position is outside of chromosome bounds");
    }

    /* set alleles at the appropriate position on the chromosome */
    ref_alleles[pos-1] = ref_nuc_id;
    alt_alleles[pos-1] = alt_nuc_id;
  }

  gzclose(gzf);
}


gzFile get_out_file(char *output_dir, int kmer_size, char *chrom_name) {
  char out_filename[1024];
  gzFile gzf;

  snprintf(out_filename, sizeof(out_filename), "%s/%s_uniq_%d.wig.gz",
	   output_dir, chrom_name, kmer_size);

  /* write wiggle header */
  gzf = util_must_gzopen(out_filename, "wb");

  fprintf(stderr, "writing output to %s\n", out_filename);

  gzprintf(gzf, "fixedStep chrom=%s start=1 step=1\n", chrom_name);

  return gzf;
}



int main(int argc, char **argv) {
  Seq *seq;
  char **chrom_names, *output_dir;
  char **fasta_files;
  GHashTable *count_tab;
  int kmer_size, n_chr, i;
  unsigned char *ref_alleles, *alt_alleles;
  gzFile gzf;
    
  if(argc < 3) {
    fprintf(stderr, "usage: %s <kmer_len> <output_dir> "
	    "<chr1.fa.gz> [<chr2.fa.gz> ...]]\n",
	    argv[0]);
    exit(2);
  }

  kmer_size = util_parse_long(argv[1]);
  output_dir = argv[2];
  n_chr = argc - 3;
  fasta_files = &argv[3];

  chrom_names = my_new(char *, n_chr);
  for(i = 0; i < n_chr; i++) {
    chrom_names[i] = get_chrom_name(fasta_files[i]);
  }
  
  seq = seq_new();

  count_tab = g_hash_table_new_full(g_str_hash, g_str_equal,
				    free, free);

  fprintf(stderr, "building count table\n");

  for(i = 0; i < n_chr; i++) {
    read_seq(seq, fasta_files[i]);

    ref_alleles = my_new(unsigned char, seq->len);
    alt_alleles = my_new(unsigned char, seq->len);
    read_snps(chrom_names[i], ref_alleles, alt_alleles, seq->len);

    fprintf(stderr, "counting kmers\n");
    count_kmers(count_tab, seq, ref_alleles, alt_alleles, kmer_size);
    
    fprintf(stderr, "counting reverse complement kmers\n");
    /* reverse complement sequence and arrays of ref and alt alleles */
    seq_revcomp(seq);
    nuc_ids_revcomp(ref_alleles, seq->len);
    nuc_ids_revcomp(alt_alleles, seq->len);
    count_kmers(count_tab, seq, ref_alleles, alt_alleles, kmer_size);

    my_free(ref_alleles);
    my_free(alt_alleles);
  }
  fprintf(stderr, "\n");  

  /* go through genome again reporting number of possible locations
   * that read maps to
   */
  fprintf(stderr, "reporting kmer counts\n");

  for(i = 0; i < n_chr; i++) {
    /* re-read sequence and SNPs for this chromosome */
    read_seq(seq, fasta_files[i]);
    ref_alleles = my_new(unsigned char, seq->len);
    alt_alleles = my_new(unsigned char, seq->len);
    read_snps(chrom_names[i], ref_alleles, alt_alleles, seq->len);

    gzf = get_out_file(output_dir, kmer_size, chrom_names[i]);
    report_kmer_hits(gzf, count_tab, seq, ref_alleles, alt_alleles,
		     kmer_size);
    gzclose(gzf);

    my_free(ref_alleles);
    my_free(alt_alleles);
  }  
  fprintf(stderr, "\n");

  fprintf(stderr, "freeing memory\n");
  for(i = 0; i < n_chr; i++) {
    my_free(chrom_names[i]);
  }
  my_free(chrom_names);

  g_hash_table_destroy(count_tab);
  seq_free(seq);
}


