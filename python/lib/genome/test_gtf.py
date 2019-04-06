import sys
import os

import chrom
import gtf


def write_gtf(filename):
    f = open(filename, "wt")

    f.write("""1	protein_coding	gene	69091	70008	.	+	.	gene_id "ENSG00000186092"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
1	protein_coding	transcript	69091	70008	.	+	.	gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547";
1	protein_coding	exon	69091	70008	.	+	.	gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547"; exon_id "ENSE00002319515";
1	protein_coding	CDS	69091	70005	.	+	0	gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547"; protein_id "ENSP00000334393";
1	protein_coding	start_codon	69091	69093	.	+	0	gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547";
1	protein_coding	stop_codon	70006	70008	.	+	0	gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS30547";
""")
    f.close()


def write_chrom_info(filename):
    f = open(filename, "wt")
    f.write("""chr1	249250621	/gbdb/hg19/hg19.2bit
chr2	243199373	/gbdb/hg19/hg19.2bit
chr3	198022430	/gbdb/hg19/hg19.2bit
chr4	191154276	/gbdb/hg19/hg19.2bit
chr5	180915260	/gbdb/hg19/hg19.2bit
chr6	171115067	/gbdb/hg19/hg19.2bit
chr7	159138663	/gbdb/hg19/hg19.2bit
chrX	155270560	/gbdb/hg19/hg19.2bit
chr8	146364022	/gbdb/hg19/hg19.2bit
chr9	141213431	/gbdb/hg19/hg19.2bit
chr10	135534747	/gbdb/hg19/hg19.2bit
chr11	135006516	/gbdb/hg19/hg19.2bit
chr12	133851895	/gbdb/hg19/hg19.2bit
chr13	115169878	/gbdb/hg19/hg19.2bit
chr14	107349540	/gbdb/hg19/hg19.2bit
chr15	102531392	/gbdb/hg19/hg19.2bit
chr16	90354753	/gbdb/hg19/hg19.2bit
chr17	81195210	/gbdb/hg19/hg19.2bit
chr18	78077248	/gbdb/hg19/hg19.2bit
chr20	63025520	/gbdb/hg19/hg19.2bit
chrY	59373566	/gbdb/hg19/hg19.2bit
chr19	59128983	/gbdb/hg19/hg19.2bit
chr22	51304566	/gbdb/hg19/hg19.2bit
chr21	48129895	/gbdb/hg19/hg19.2bit
""")
    f.close()            

            
            
class TestGTF:

    def get_chrom_dict(self):
        filename = "test_chrom_info.txt"
        write_chrom_info(filename)
        chrom_dict = chrom.parse_chromosomes_dict(filename)
        os.unlink(filename)
        return chrom_dict


    def test_gtf(self):
        gtf_filename = "test_gtf.txt"
        write_gtf(gtf_filename)

        chrom_dict = self.get_chrom_dict()
        gene_list, gene_dict, tr_list, tr_dict = gtf.parse_gtf(gtf_filename, chrom_dict)

        assert(len(gene_list) == 1)

        g = gene_dict['ENSG00000186092']

        assert(g.start == 69091)
        assert(g.end == 70008)

        assert(g.id == "ENSG00000186092")
        assert(g.name == "OR4F5")

        assert(g.biotype == "protein_coding")
        assert(g.source == "ensembl_havana")

        assert(len(g.transcripts) == 1)

        tr = g.transcripts[0]
        assert(tr.cds_start == 69091)
        assert(tr.cds_end == 70005)

        assert(len(tr.exons) == 1)
        assert(tr.exons[0].start == 69091)
        assert(tr.exons[0].end == 70008)

        tr = g.transcripts[0]

        ### TODO: check transcript and exons
        
        os.unlink(gtf_filename)
