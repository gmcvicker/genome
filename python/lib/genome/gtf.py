import sys

import genome.coord
import genome.gene
import genome.transcript


def parse_gtf(path, chrom_dict):
    """Reads a list of genes from a GTF file"""
    f = open(path, "rt")
    
    gene_dict = {}
    gene_list = []
    tr_dict = {}

    for line in f:
        if line.startswith("#"):
            # comment line
            continue

        words = line.split("\t")

        feature_type = words[2]
        chrom = parse_chrom(words[0], chrom_dict)
        
        if feature_type == "gene":
            g = parse_gene(words, chrom)

            if g.id in gene_dict:
                raise ValueError("multiple genes with ID '%s'" % g.id)
            gene_dict[g.id] = g
            gene_list.append(g)
            

        elif feature_type == "transcript":
            tr = parse_transcript(words, chrom)
            if tr.id in tr_dict:
                raise ValueError("Multiple transcripts with ID '%s'" %
                                 tr.id)
            tr_dict[tr.id] = tr

            if tr.gene_id in gene_dict:
                gene_dict[tr.gene_id].add_transcript(tr)
            else:
                # for now assume that gene will come before
                # transcript in GTF file, not sure this is
                # guaranteed to be true...
                raise ValueError("expected gene with id %s to come before "
                                 "transcript %s in GTF file" %
                                 (tr.gene_id, tr.id))
                    
        elif feature_type == "exon":
            exon = parse_exon(words, chrom)

            if exon.transcript_id in tr_dict:
                tr_dict[exon.transcript_id].add_exon(exon)
            else:
                # for now assume that transcript will come before
                # exon in GTF file, not sure this guaranteed to be true...
                raise ValueError("expected gene with id %s to come before "
                                 "transcript %s in GTF file" %
                                 (tr.gene_id, tr.id))
                
        elif feature_type == "CDS":
            parse_cds(words, tr_dict)
        else:
            # do nothing with other annotation types
            pass

    # also return transcript list and transcript dict????
        
    return gene_list, gene_dict
            
                  
def parse_attrib(attr_str):
    """parse ';'-delimited string from attribute column of a GTF file"""
    attr_dict = {}
    attr_str = attr_str.strip()
    if attr_str[-1] == ";":
        # remove trailing ';'
        attr_str = attr_str[:-1]
    
    for attr in attr_str.strip().split(";"):
        attr = attr.strip()
        try:
            key, val = attr.split(" ")
        except:
            raise ValueError("failed to parse attribute string: '%s'\n",
                             attr)
        
        val = val.replace('"', '')
        attr_dict[key] = val    
    return attr_dict


def parse_chrom(chrom_str, chrom_dict):
    # try adding "chr" to start of chromosome
    # if the name is not found
    if chrom_str not in chrom_dict:
        new_name = "chr" + chrom_str
        if new_name in chrom_dict:
            chrom_str = new_name
        else:
            raise genome.coord.CoordError("unknown chromosome '%s'" %
                                          chrom_str)
    return chrom_dict[chrom_str]
    


def parse_transcript(words, chrom):
    start = int(words[3])
    end = int(words[4])
    strand = genome.coord.parse_strand(words[6])

    attr_dict = parse_attrib(words[8])

    gene_id = None
    if 'gene_id' in attr_dict:
        gene_id = attr_dict['gene_id']

    tr_id = None
    if 'transcript_id' in attr_dict:
        tr_id = attr_dict['transcript_id']

    tr_name = None
    if 'transcript_name' in attr_dict:
        tr_name = attr_dict['transcript_name']
        
        
    tr = genome.transcript.Transcript(chrom=chrom, name=tr_name, start=start,
                                      end=end, strand=strand, id=tr_id,
                                      gene_id=gene_id)

    return tr
    
    

def parse_gene(words, chrom):
    gene_type = words[1]
    start = int(words[3])
    end = int(words[4])
    score = None if(words[5] == '.') else float(words[5])

    strand = genome.coord.parse_strand(words[6])

    attr_dict = parse_attrib(words[8])
    
    id_str = attr_dict['gene_id'] if 'gene_id' in attr_dict else None

    if 'gene_biotype' in attr_dict:
        biotype = attr_dict['gene_biotype']
    else:
        biotype = None

    if 'gene_name' in attr_dict:
        name = attr_dict['gene_name']
    else:
        name = None

    if 'gene_source' in attr_dict:
        source = attr_dict['gene_source']
    else:
        source = None

    g = genome.gene.Gene(chrom=chrom, start=start, end=end, strand=strand,
                         score=score, id=id_str, biotype=biotype, source=source)

    return g
    



def parse_exon(words, chrom):
    start = int(words[3])
    end = int(words[4])
    strand = genome.coord.parse_strand(words[6])

    attr_dict = parse_attrib(words[8])

    tr_id = None
    if 'transcript_id' in attr_dict:
        tr_id = attr_dict['transcript_id']

    exon = genome.coord.Coord(chrom, start, end, strand=strand)
    exon.transcript_id = tr_id

    return exon
    
    
def parse_cds(words, tr_dict):
    start = int(words[3])
    end = int(words[4])

    attr_dict = parse_attrib(words[8])

    tr_id = None
    if 'transcript_id' in attr_dict:
        tr_id = attr_dict['transcript_id']

    if tr_id in tr_dict:
        tr = tr_dict[tr_id]

        # update CDS start/end
        if (tr.cds_start is None) or (tr.cds_start > start):
            tr.cds_start = start
        if (tr.cds_end is None) or (tr.cds_end < end):
            tr.cds_end = end
    else:
        raise ValueError("expected gene with id %s to come before "
                         "transcript %s in GTF file" %
                         (tr.gene_id, tr.id))

            
        
        
    
