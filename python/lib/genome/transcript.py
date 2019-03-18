#!/usr/bin/env python

import sys

from genome.coord import Coord, CoordError

class Transcript(Coord):
    def __init__(self, chrom=None, name=None, exons=[],
                 start=None, end=None, strand=None, id=None, gene_id=None,
                 cds_start=None, cds_end=None):
        """Creates and initializes a new transcript object."""

        self.chrom = chrom # chromsome object
        self.start = start # start of this transcript
        self.end = end # end of this transcript
        self.strand = strand # strand of this transcript
        self.id = id # identifier (e.g. ENST identifier) for this transcript
        self.gene_id = gene_id # id of gene this transcript associated with
        self.exons = [] # exons for this transcript

        for ex in exons:
            self.add_exon(ex)

        self.name = name

        # genomic coordinates of coding sequence start/end
        # set to None if non-coding transcript
        self.cds_start = cds_start
        self.cds_end = cds_end

        # self.check_exon_coords()

    
    def n_exon(self):
        """Returns the number of exons in this transcript"""
        return len(self.exons)

    
    def size(self):
        """Gives the combined length of exons in this transcript
        (cDNA length)"""
        ttl_size = 0
        for exon in self.exons:
            ttl_size += exon.length()
        
        return ttl_size


    def is_coding(self):
        """Returns true if this transcript is protein-coding, FALSE
        otherwise"""
        return self.cds_start is not None


    def copy(self):
        """Returns a copy of this transcript. The exons are all copied,
        but the chromosomes pointed to by the exons are not."""
        new_exons = []
        for exon in self.exons:
            new_exons.append(exon.copy())

        return Transcript(chrom=self.chrom, name=self.name, exons=new_exons,
                          start=self.start, end=self.end, strand=self.strand,
                          id=self.id, gene_id=self.gene_id,
                          cds_start=self.cds_start,
                          cds_end=self.cds_end)
                          



    def __check_transcript_extent(self):
        """checks that transcript coordinates appear consistant with
        region spanned by first/last exons"""
        if len(self.exons) == 0:
            return
        
        first_exon = self.exons[0]
        last_exon = self.exons[-1]

        if first_exon.strand != self.strand:
            raise CoordError("transcript strand is not consistant"
                                  "with exon strand")

        if self.strand == 1:
            if first_exon.start != self.start or \
               last_exon.end != self.end:
                raise CoordError("transcript coordinates are not "
                                      "consistant with first/last exon"
                                      "coordinates")
        elif self.strand == -1:
            if first_exon.end != self.end or \
               last_exon.start != self.start:
                raise CoordError("transcript coordinates (%s) are not "
                                      "consistant with first (%s) and "
                                      "last (%s) exon coordinates" % \
                                      (str(self), str(first_exon),
                                       str(last_exon)))


    def __check_exon_ordering(self):
        """checks that exon coordinates are consistant with each other in
        ordering, chromosome, and strand"""
        if len(self.exons) == 0:
            return
        
        prev_exon = self.exons[0]
        exon_num = 1
        
        for exon in self.exons[1:]:
            exon_num += 1
            if exon.strand != prev_exon.strand:
                raise CoordError("exon strands are inconsistant")
            
            if exon.chrom.idnum != prev_exon.chrom.idnum:
                raise CoordError("exon chromosomes are inconsistant")
            
            if self.strand == 1:
                # each exon must be separated from
                # previous exon by intron of at least 1 bp
                if exon.start <= prev_exon.end+1:
                    raise CoordError("exon ordering or coordinates "
                                          "are not consistant:\n" +
                                          "  exon %d: %s\n  exon %d: %s" % \
                                          (exon_num-1, str(prev_exon),
                                           exon_num, str(exon)))
            else:
                if exon.end >= prev_exon.start-1:
                    raise CoordError("exon ordering or coordinates "
                                          "are not consistant:\n"                                                            "  exon %d:%s\n  exon %d:%s" % \
                                          (exon_num-1, str(prev_exon),
                                           exon_num, str(exon)))


    def __check_cds_coords(self):
        """verifies that the cds start and end of this transcript
        fall into exons and are ordered correctly"""

        if self.cds_start is None and self.cds_end is not None:
            raise CoordError("Transcript defines CDS start but not end")
        if self.cds_end is None and self.cds_start is not None:
            raise CoordError("Transcript defines CDS end but not start")
        
        if (self.cds_start is None and self.cds_end is None):
            return
        
        # CDS must span at least 3 bp
        span = self.cds_end - self.cds_start
        if span < 3:
            raise CoordError("Expected region spanned by"
                                  "CDS to be at least 3 bp, got %d" % span)
        
        # check that CDS start and end positions are contained within exons
        cds_start_in_exon = False
        cds_end_in_exon = False

        for ex in self.exons:
            if self.cds_start >= ex.start and self.cds_start <= ex.end:
                cds_start_in_exon = True
            if self.cds_end >= ex.start and self.cds_end <= ex.end:
                cds_end_in_exon = True

        if not cds_start_in_exon:
            raise CoordError("CDS start position not in exon")
        if not cds_end_in_exon:
            raise CoordError("CDS end position not in exon")

    
    def check_exon_coords(self):
        """Checks the coordinates of the exons of this
        transcript to make sure they are ordered correctly, etc.
        Raises a CoordError if a problem is detected"""
        if self.n_exon() < 1:
            return

        self.__check_transcript_extent()
        self.__check_exon_ordering()
        self.__check_cds_coords()

    
    def add_exon(self, exon):
        """Adds an exon to this transcript. If the exon falls outside of the
        coordinates of this transcript, the transcript
        coordinates are updated"""
        
        if self.n_exon() == 0:
            # this is first exon added to transcript
            self.chrom = exon.chrom
            self.start = exon.start
            self.end = exon.end
            self.strand = exon.strand
        else:
            if exon.chrom.idnum != self.chrom.idnum:
                raise ValueError("exon chromosome %s does not match "
                                 "transcript chromosome %s" % \
                                 (str(exon.chrom), str(self.chrom)))
            
            if exon.strand != self.strand:
                raise ValueError("exon strand %d does not match "
                                 "transcript strand %d" % (exon.strand,
                                                           self.strand))
            
            if exon.start < self.start:
                self.start = exon.start
            if exon.end > self.end:
                self.end = exon.end
                
        self.exons.append(exon)


    def get_introns(self):
        """Returns a list of coordinates that represent the introns for
        this transcript"""
        introns = []

        for i in range(len(self.exons)-1):
            ex1 = self.exons[i]
            ex2 = self.exons[i+1]

            if self.strand == -1:
                intron = Coord(self.chrom, ex2.end+1, ex1.start-1,
                               strand=self.strand)
            else:
                intron = Coord(self.chrom, ex1.end+1, ex2.start-1,
                               strand=self.strand)

            intron.exon_5p = ex1
            intron.exon_3p = ex2
            introns.append(intron)
            
        return introns



    def __str__(self):
        id_str = "NA"
        if self.id is not None:
            id_str = self.id
        
        name_str = "NA"
        if self.name is not None:
            name_str = self.name

        strand_str = "0"
        if self.strand is not None:
            strand_str = str(self.strand)

        cds_start_str = "NA"
        if self.cds_start is not None:
            cds_start_str = str(self.cds_start)

        cds_end_str = "NA"
        if self.cds_end is not None:
            cds_end_str = str(self.cds_end)

        exon_starts = [str(ex.start) for ex in self.exons]
        exon_ends = [str(ex.end) for ex in self.exons]

        exon_start_str = ",".join(exon_starts)
        exon_end_str = ",".join(exon_ends)

        fields = [id_str, name_str, self.chrom.name,
                  str(self.start), str(self.end), strand_str,
                  exon_start_str, exon_end_str, cds_start_str,
                  cds_end_str]

        return "\t".join(fields)



    def update_bounds(self):
        """Updates the start and end of this transcript to reflect the
        start / end of the first / last exons. This is useful if the
        exons have been updated."""
        if self.strand == -1:
            left_exon = self.exons[-1]
            right_exon = self.exons[0]
        else:
            left_exon = self.exons[0]
            right_exon = self.exons[-1]

        if left_exon.start != self.start:
            self.start = left_exon.start
        if right_exon.end != self.end:
            self.end = right_exon.end
        
#
# TODO:
# need to update read_transcripts so that it can read from GTF file
def read_transcripts(path, chrom_dict):
    """Retrives all transcripts from the specified transcript file"""

    f = open(path, "r")

    transcripts = []

    header = f.readline().rstrip().split()
    row = dict([(header[x], x) for x in range(len(header))])
    
    for line in f:
        row = dict(zip(header, line.rstrip().split()))
        if row['ID'] == "NA":
            tr_id = None
        else:
            tr_id = row['ID']

        if row["NAME"] == "NA":
            name = None
        else:
            name = row['NAME']
        
        # parse CDS start/end
        if row['CDS.START'] == 'NA':
            cds_start = None
        else:
            cds_start = int(row['CDS.START'])

        if row['CDS.END'] == 'NA':
            cds_end = None
        else:
            cds_end = int(row['CDS.END'])

        strand = int(row['STRAND'])
        chrom = chrom_dict[row['CHROM']]

        # parse exons
        exon_starts = [int(x) for x in row['EXON.STARTS'].split(",")]
        exon_ends = [int(x) for x in row['EXON.ENDS'].split(",")]

        exons = []
        for i in range(len(exon_starts)):
            exon = Coord(chrom, exon_starts[i], exon_ends[i], strand)
            exons.append(exon)
        
        tr = Transcript(name=name, exons=exons,
                        cds_start=cds_start, cds_end=cds_end,
                        id=tr_id)

        transcripts.append(tr)
        
    f.close()

    return transcripts




if __name__ == "__main__":
    pass
    # if len(sys.argv) != 2:
    #     sys.stderr.write("usage: %s <transcript_file>\n" % sys.argv[0])
    #     exit(255)

    # with genome.db.connect(mode="r+") as gdb:
    #     chrom_adp = gdb.get_adaptor("chromosome")
    #     chrom_dict = chrom_adp.fetch_name_dict()
 
    # trs = read_transcripts(sys.argv[1], chrom_dict)

    # for tr in trs:
    #     print((str(tr)))

