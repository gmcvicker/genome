import sys

from genome import coord
import numpy as np

class Gene(coord.Coord):
    def __init__(self, transcripts=[], score=None, idnum=None):
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.idnum = idnum
        self.score = score
        self.transcripts = []

        try:
            for tr in transcripts:
                self.add_transcript(tr)
        except ValueError:
            for tr in transcripts:
                sys.stderr.write("  " + str(tr) + "\n")
                raise ValueError("invalid transcripts for gene")

    def add_transcript(self, tr):
        """Adds a transcript"""
        if self.chrom is None:
            self.chrom = tr.chrom
            self.start = tr.start
            self.end = tr.end
            self.strand = tr.strand
        else:
            # check transcript is on same strand and chromosome
            if self.chrom.name != tr.chrom.name:
                raise ValueError("mismatched chromosomes")
            if self.strand != tr.strand:
                raise ValueError("mismatched strands")

            # update start/end
            if self.start > tr.start:
                self.start = tr.start
                
            if self.end < tr.end:
                self.end = tr.end

        self.transcripts.append(tr)


    def get_merged_exons(self):
        """Returns a list of non-redundant exon coordinates for this
        gene. Overlapping exons are merged."""

        # get redundant, sorted list of all exons
        exon_list = []
        for tr in self.transcripts:
            exon_list.extend(tr.exons)
        coord.sort_coords(exon_list, use_strand=False)

        # create coordinate group from first exon
        cur_exon = coord.CoordGroup(exon_list[0])
        merged_exons = [cur_exon]

        # merge overlapping exons
        for ex in exon_list[1:]:
            if ex.overlaps(cur_exon):
                # expand current exon
                cur_exon.add_coord(ex)
            else:
                # new exon group
                cur_exon = coord.CoordGroup(ex)
                merged_exons.append(cur_exon)

        return merged_exons


    def get_longest_transcript(self):
        max_sz = 0
        max_tr = None
        for tr in self.transcripts:
            sz = tr.size()
            if sz > max_sz:
                max_tr = tr
                max_sz = sz

        return max_tr
                

    def get_unique_tss(self):
        """Returns a list of unique transcription start site coordinates
        for this gene"""
        tss_dict = {}

        for tr in self.transcripts:
            if tr.strand == 1:
                if tr.start in tss_dict:
                    pass
                else:
                    tss = coord.Coord(self.chrom, tr.start, tr.start,
                                      strand=self.strand)
                    tss_dict[tss.start] = tss
            else:
                if tr.end in tss_dict:
                    pass
                else:
                    tss = coord.Coord(self.chrom, tr.end, tr.end,
                                      strand=self.strand)
                    tss_dict[tss.end] = tss

        return tss_dict.values()
                    


    def get_unique_exons(self):
        """Returns list of exons with unique start/end coordinates
        for this gene, discarding duplicates with the same start, end,
        and strand"""
        exon_list = []
        seen_exons = set([])
        for tr in self.transcripts:
            for ex in tr.exons:
                key = ":".join([str(ex.start), str(ex.end), str(ex.strand)])
                if key not in seen_exons:
                    seen_exons.add(key)
                    exon_list.append(ex)
        return exon_list



    def get_unique_introns(self):
        """Returns list of non-redundant introns from all of the
        transcripts in this gene"""
        seen_introns = set([])
        intron_list = []
        for tr in self.transcripts:
            for intron in tr.get_introns():
                key = ":".join([intron.chrom.name, str(intron.start),
                                str(intron.end)])
                if key in seen_introns:
                    pass
                else:
                    seen_introns.add(key)
                    intron_list.append(intron)
        return intron_list



    def mask_exons(self, mask, mask_region):
        """Given a mask, and a coordinate that defines the region
        spanned by the mask, sets all portions of the mask spanned by
        this gene's exons to True"""
        for ex in self.get_merged_exons():
            # convert exon coords into mask slice:
            # 0 is the start of the mask region
            start = ex.start - mask_region.start
            end = ex.end - mask_region.end + 1
            mask[start:end] = True

    def get_exon_mask(self):
        """Returns a boolean numpy array the length of the gene
        with every exonic position set to True, and every non-exonic
        position set to False."""
        mask = np.zeros(self.length(), dtype=np.bool)
        self.mask_exons(mask, self)
        return mask


    def update_bounds(self):
        """Updates the start/end of all transcripts associated with this
        gene and then updates the start/end of this gene. This is useful
        if the start/end positions of the exons have changed."""
        min_start = None
        max_end = None
        for tr in self.transcripts:
            tr.update_bounds()
            if min_start is None or tr.start < min_start:
                min_start = tr.start
            if max_end is None or tr.end > max_end:
                max_end = tr.end

            self.start = min_start
            self.end = max_end
            




def group_transcripts(trs):
    """Creates a list of genes created from overlapping sets of
    transcripts. The provided list of transcripts is sorted in-place
    by this function."""

    if len(trs) == 0:
        return []
    
    # sort transcripts, and then find overlapping transcripts on same strand
    coord.sort_coords(trs, use_strand=True)

    # start first gene with first transcript
    cur_gene = Gene([trs[0]])
    genes = [cur_gene]

    for tr in trs[1:]:
        if tr.overlaps(cur_gene, use_strand=True):
            # keep adding overlapping transcripts to current gene
            cur_gene.add_transcript(tr)
        else:
            # start a new gene
            cur_gene = Gene([tr])
            genes.append(cur_gene)

    return genes



