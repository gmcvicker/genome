#!/usr/bin/env python

import coord

class Alignment(coord.Coord):
    def __init__(self, crossmatch_header=None, query_seq=None, target_seq=None):
        self.query_seq = query_seq
        self.target_seq = target_seq
        self.coverage = None
        
        if not crossmatch_header is None:
            self.parse_header(crossmatch_header)
        
            self.ID = self.query_id
        
            # make a unique read id
            self.read_id = self.chr + ":" + self.strand + ":" + \
                           str(self.start) + ":" + str(self.end) + ":" + str(self.query_start) + \
                           ":" + str(self.query_end)
    
    def parse_header(self, header):
        words = header.rstrip().rsplit()
        self.score = int(words[1])
        self.percent_diff = float(words[2])
        self.percent_query_gap = float(words[3])
        self.percent_target_gap = float(words[4])
        self.query_id = words[5]
        self.query_start = int(words[6])
        self.query_end = int(words[7])
        self.query_past_end = words[8]
        
        if words[9] == "C":
            self.chr = words[10]
            self.strand = "-"
            self.start = int(words[13])
            self.end = int(words[12])
            self.chr_past_end = words[11]
        else:
            self.chr = words[9]
            self.strand = "+"
            self.start = int(words[10])
            self.end = int(words[11])
            self.chr_past_end = words[12]
            
        if self.start > self.end:
            raise Exception("chromosome start should not be greater than end")
            
        if self.query_start > self.query_end:
            raise Exception("query start should not be greater than end")
            
    def __str__(self):
        if self.strand == "+":            
            coord_strs = [str(self.query_start), str(self.query_end), self.query_past_end,
                          self.chr, str(self.start), str(self.end), str(self.chr_past_end)]
        else:
            coord_strs = [str(self.query_start), str(self.query_end), self.query_past_end,
                          "C", self.chr, str(self.chr_past_end), str(self.end), str(self.start)]
        
        align_str = ""
        if self.query_seq and self.target_seq:
            # only include alignment lines if they are defined
            align_str = "\n" + self.query_seq + "\n" + self.target_seq + "\n"
        
        return "\t".join(["ALIGNMENT", str(self.score), str(self.percent_diff),
                          str(self.percent_query_gap), str(self.percent_target_gap),
                          self.query_id] + coord_strs) + align_str

    def __cmp__(self, other):
        """A comparison function, used by the sort function to decide
        the order of Coords. If the strand attribute is not None, it
        is used as part of the comparison."""
        
        return cmp(self.chr, other.chr) or cmp(self.start, other.start)


def parse_crossmatch_seqs(file, alignment):
    # get the first line of alignments
    line = file.readline()
    
    query_seq = ""
    target_seq = ""
    
    while line:
        # stop once we hit "Transitions" line
        if line.startswith("Transitions"):
            break
        
        # current line should be query line
        words = line.rstrip().rsplit()
        if words[0] == "C":
            # append to query sequence
            query_seq += words[3]
        else:
            query_seq += words[2]
        
        # next line is difference line, we don't need it
        line = file.readline()
        # need to read again?? there is some weird \r or something
        line = file.readline()
        
        # next line contains target sequence
        words = line.rstrip().rsplit()
        target_seq += words[2]
        
        # skip next blank line and read next query line
        line = file.readline()
        line = file.readline()
        
    alignment.query_seq = query_seq
    alignment.target_seq = target_seq



def parse_seqs(file, alignment):
    # next line is a blank line
    line = file.readline().rstrip()
    
    if line.startswith("ALIGNMENT"):
        raise exception("No sequences found in alignment")
    
    if line == "":
        # blank line indicates that
        # this is a crossmatch-style alignment where seqs are
        # split across multiple lines
        parse_crossmatch_seqs(file, alignment)
    else:
        # this is a simple alignment where seqs are one after each other
        alignment.query_seq = line
        alignment.target_seq = file.readline().rstrip()


def parse_crossmatch(file, include_seqs=True):
    """Parse an alignment block that starts with an ALIGNMENT line and ends with
    a Transitions / transversions line"""
    
    # advance to next ALIGNMENT line:
    line = file.readline()
    while line:
        if line.startswith("ALIGN"):
            header_line = line.rstrip()
            break
        line = file.readline()
    
    if line == "":
        # we are at end of file
        return None
    
    # create alignment from the header line
    alignment = Alignment(header_line)

    if include_seqs:
        parse_seqs(file, alignment)
    
    return alignment

