
from chrom import Chromosome
import math
import sys
import gzip

import numpy as np

class CoordError(Exception):
    """An exception indicating that something is wrong with a coordinate"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(self.value)


class Coord(object):
    def __init__(self, chrom, start, end, strand=0,
                 score=None, idnum=None, name=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.idnum = idnum
        self.score = score
        self.name = name

        if start > end:
            raise CoordError("start (%d) should be less than or "
                                  "equal to end (%d)" % (start, end))

        if start < 1:
            raise CoordError("start (%d) should not be less "
                                  "than 1" % start)

        if end > chrom.length:
            raise CoordError("end (%d) should not be greater than "
                                  "length of chromosome "
                                  "(%d)" % (end, chrom.length))
        
        if strand != 0 and strand != 1 and strand != -1:
            raise CoordError("strand should be one of (-1, 0, 1)")


    def key(self, use_strand=False):
        """Returns a tuple to be used for sorting of coordinates. If
        use_strand is True the tuple consists of (chrom.idnum, strand, start),
        otherwise the tuple is (chrom.idnum, start)"""
        if use_strand:
            return (self.chrom.idnum, self.strand, self.start)

        return (self.chrom.idnum, self.start)
        
        
    def __str__(self):
        """returns a string representation of this coordinate"""
        if self.idnum:
            id_str = str(self.idnum) + " "
        else:
            id_str = ""
            
        if self.strand == 1:
            strand_str = "(+)"
        elif self.strand == -1:
            strand_str = "(-)"
        else:
            strand_str = "(.)"
        
        return(id_str + str(self.chrom) + ":" + str(self.start) + "-" + \
               str(self.end) + strand_str)
    
    
    def length(self):
        """Returns the size of the region spanned by the
        coordinate in bases"""
        return self.end - self.start + 1
    
    def overlaps(self, other, use_strand=False):
        """Returns True if the coords overlap and False if they don't.
        If use_strand is True then the coordinates must be on the same
        strand to be considered overlapping."""
        if self.chrom.idnum != other.chrom.idnum:
            return False
        if use_strand and self.strand != other.strand:
            return False
        if self.start <= other.end and self.end >= other.start:
            return True
        return False

    def copy(self):
        """Creates a copy of this coordinate object and returns it.
        The copy is shallow in the sense that only the reference to
        the chromosome attribute is copied (i.e. a new Chromosome
        is not created)."""
        return Coord(self.chrom, self.start, self.end,
                     strand=self.strand, score=self.score,
                     idnum=self.idnum)


    def expand(self, n_bp):
        """Expands this coordinate by n_bp in both directions, but not
        exceeding the boundaries of the chromosome that this coordinate
        is on"""
        self.start = self.start - n_bp
        self.end = self.end + n_bp

        if self.chrom:
            # if there is an associated chromosome, don't go past the ends
            if self.start < 1:
                self.start = 1
            if self.end > self.chrom.length:
                self.end = self.chrom.length

    def within(self, other, use_strand=False):
        """checks whether one coordinate is completely contained
        within the other one. If use_strand is True then the
        this coordinate must also be on the same strand as the
        other one, in order to be considered 'contained'"""
        if self.chrom and other.chrom and (self.chrom.name != other.chrom.name):
            return False
        
        if use_strand and self.strand != other.strand:
            return False

        if self.start >= other.start and self.end <= other.end:
            return True

        return False


class CoordGroup(Coord):
    def __init__(self, coord):
        """Initializes CoordGroup object from a coordinate object"""
        # call superclass constructor:
        super(CoordGroup, self).__init__(coord.chrom, coord.start,
                                         coord.end, strand=coord.strand,
                                         score=coord.score,
                                         idnum=coord.idnum,
                                         name=coord.name)
        self.coord_list = [coord]
        
    def add_coord(self, coord, force=False, use_strand=False):
        """adds a read to this read group, expands the read group as
        necessary. If the force flag is true, the read is allowed to be
        added even if it does not overlap"""
        
        if not force and (coord.end < self.start or coord.start > self.end):
            raise CoordError("Coord does not overlap "
                                  "coordinate group")
        
        if use_strand and coord.strand != self.strand:
            raise CoordError("Coord is not on same strand as "
                                  "CoordGroup")
        
        if coord.chrom.idnum != self.chrom.idnum:
            raise CoordError("Coord is not on same chrom as "
                                  "CoordGroup")
        
        if coord.start < self.start:
            # update start coordinate of this read group
            self.start = coord.start
        if coord.end > self.end:
            # update end of coordinate of this read group
            self.end = coord.end
        
        self.coord_list.append(coord)
        
    def num_coords(self):
        """returns the number of coordinates in this coordinate group"""
        return len(self.coord_list)
    
    
    def copy(self):
        """Creates a read group which is a copy of this one"""
        new_group = CoordGroup(self)
        new_group.coord_list = list(self.coord_list)
        return new_group
    
    def split(self, max_coords_per_group):
        """Returns a list of read groups with fewer reads in each group"""
        if self.num_coords() < max_coords_per_group:
            return [self.copy()]
        
        # determine the number of groups we need to make
        n_groups = int(math.ceil(float(self.num_coords()) /
                                 float(max_coords_per_group)))
        
        coords_per_group = int(math.ceil(float(self.num_coords()) /
                                         float(n_groups)))
        
        # now create new read groups
        new_groups = []
        cur_coord_group = None
        cur_n_coord = 0
        cur_group_id = 0
        for coord in self.coord_list:
            if cur_n_coord >= coords_per_group or cur_coord_group is None:
                # create a new coordinate group with this coord
                cur_coord_group = CoordGroup(coord)
                cur_group_id += 1
                
                # rename cur group id to give the group number at end
                cur_coord_group.idnum = self.idnum + "_RG" + str(cur_group_id)
                new_groups.append(cur_coord_group)
                cur_n_coord = 1
            else:
                # add coordinate to existing coordinate group
                if cur_coord_group.overlaps(coord):
                    cur_coord_group.add_coord(coord)
                else:
                    # unusual situation in which reads do not overlap
                    # even though they are ordered by start position
                    # this can happen when we split read groups and
                    # the aligned portion of one of the reads is shorter
                    # than the full read length
                    cur_coord_group.add_coord(coord, force=True)
                    sys.stderr.write("WARNING: adding non-overlapping "
                                     "read to split read group\n")
                    
                cur_n_coord += 1
        
        return new_groups
        
        

def parse_strand(strand_str):
    """Parses a strand string that can be in several possible formats.
    Returns 1 for forward strand, -1 for reverse strand and 0 for
    unknown/undefined strand."""
    if(strand_str in ('+', '1', 'f', 'F', 'fwd', 'forward')):
        return 1

    if(strand_str in ('-', '-1', 'r', 'R', 'rev', 'reverse')):
        return -1

    if(strand_str in (".", '0', 'u', 'U', '')):
        return 0

    raise CoordError("unknown strand: '%s'" % strand_str)

        


def sort_coords(coords, use_strand=False):
    """Sorts the provided list of coordinates in-place. If use_strand
    is True, then coordinates are sorted first by chromosome then
    strand, then start position (otherwise by chromosome, and start
    position)"""

    # define a closure, which takes only a single arg
    # but will remember whether or not to use strand
    def key_func(coord):
        return Coord.key(coord, use_strand=use_strand)

    coords.sort(key=key_func)


def get_coord_overlaps(coord, coord_list, use_strand=False):
    """Wrapper for get_overlaps function that given a single coordinate
    and a list of sorted coordinates, returns a list of all coordinates
    from the list that overlap with single coordinate"""
    overlap_list = get_overlaps([coord], coord_list, use_strand=False)
    return overlap_list[0]




def coords_from_sites(chrom, sites):
    """Creates a list of Coord objects from a list of sites.
    Sites that are adjacent are combined into a single coordinate.
    For example the position list [1,2,10,11,12,13,20] would give
    a list of coordinates spanning the following ranges
    [1-2, 10-13, 20-20]."""
    start = None
    end = None

    coords = []

    for site in sites:
        if start is None:
            # start first block
            end = start = site
        else:
            if (site - end) != 1:
                # write out previous block
                coords.append(Coord(chrom, start, end))
                # start new block
                start = site
            # record end of current block
            end = site

    if end:
        # finish last block
        coords.append(Coord(chrom, start+1, end+1))

    return coords
    


def get_overlaps(coords1, coords2, use_strand=False):
    """Takes two sorted lists of coordinates and returns a list with indices
    corresponding to the first list of coordinates. Any coordinates from
    the second list, which overlap the first are provided as lists at each
    index. For example the following represents overlaps of coordinates i,
    i+1, and i+2:
      [i]   => [coord2_x, coord2_y]
      [i+1] => [coord2_y]
      [i+2] => []
    
    If use_strand is True then only coordinates that are on the same
    strand are considered overlapping (note: the provided coordinates
    must have been sorted using strand information in this case)."""
    i = 0
    j = 0

    # for brevity:
    s = use_strand

    # create list of empty lists
    overlap_list = [[] for x in xrange(0, len(coords1))]
    
    # loop over both lists until the end of one of them reached
    while i < len(coords1) and j < len(coords2):
        # move through first list until you overlap or pass the current
        # coordinate in the second list
        while (i < len(coords1)) and \
            (not coords1[i].overlaps(coords2[j], use_strand=use_strand)) and \
            (cmp(coords1[i].key(use_strand=s),
                 coords2[j].key(use_strand=s)) < 0):
            i += 1
        
        # move through second list until you overlap or pass the current
        # coordinate in the first list
        while (i < len(coords1)) and (j < len(coords2)) and \
            (not coords2[j].overlaps(coords1[i], use_strand=s)) and \
            (cmp(coords2[j].key(use_strand=s),
                 coords1[i].key(use_strand=s)) < 0):
            j += 1
        
        # keep adding overlapping coords, but keep track of where we
        # were in read list, because same coord can overlap with many
        j_overlap = j        
        while (j_overlap < len(coords2)) and (i < len(coords1)) and \
            (cmp(coords2[j_overlap].key(use_strand=s),
                 coords1[i].key(use_strand=s)) < 0 or \
             coords2[j_overlap].overlaps(coords1[i], use_strand=s)):

            if coords2[j_overlap].overlaps(coords1[i], use_strand=s):
                overlap_list[i].append(coords2[j_overlap])

            j_overlap += 1

        i += 1
        
    return overlap_list







def np_overlap_strand(c1, c2):
    """Returns True if the provided numpy coordinate array elements overlap
    and are on the same strand, returns False otherwise"""
    return((c1['chromosome_id'] == c2['chromosome_id']) and
           (c1['start']  <= c2['end']) and
           (c1['end']    >= c2['start']) and
           (c1['strand'] == c2['strand']))


def np_overlap(c1, c2):
    """Returns True if the numpy coordinate array elements overlap
    (ignoring strand), returns False otherwise"""
    return((c1['chromosome_id'] == c2['chromosome_id']) and
           (c1['start'] <= c2['end']) and
           (c1['end']   >= c2['start']))


def np_cmp(c1, c2):
    """Performs a comparison of two numpy coordinate array elements.
    Returns -1 if c1 is less than c2, 0 if they are the same, and
    1 if c1 is greater than c2. The comparison is performed numerically
    and uses first the chromosome_id and then the start"""
    return cmp((c1['chromosome_id'], c1['start']),
               (c2['chromosome_id'], c2['start']))


def np_cmp_strand(c1, c2):
    """Performs a comparison of two numpy coordinate array elements,
    as np_cmp function, but including strand as the second element
    of the comparison: (chromosome_id, strand, start)"""
    return cmp((c1['chromosome_id'], c1['strand'], c1['start']),
               (c2['chromosome_id'], c2['strand'], c2['start']))
        


def get_np_overlaps(coords1, coords2, use_strand=False):
    """Takes two sorted numpy coordinate arrays and returns a list with indices
    corresponding to the first array of coordinates. Coordinates from
    the second list, which overlap the first are provided as numpy arrays
    at each index. For example the following represents overlaps of
    coordinates i, i+1, and i+2:
      [i]   => [coord2_x, coord2_y]
      [i+1] => [coord2_y]
      [i+2] => []    
    If use_strand is True then only coordinates that are on the same
    strand are considered overlapping (note: the provided coordinates
    must have been sorted using (chromosome_id, strand, start) as an
    ordering if use_strand is True))."""
    i = 0
    j = 0    

    if use_strand:
        ovlp = np_overlap_strand
        compare = np_cmp_strand
    else:
        ovlp = np_overlap
        compare = np_cmp

    # create list of empty lists
    overlap_list = [[] for x in xrange(0, coords1.size)]
    
    # loop over both lists until the end of one of them reached
    while i < coords1.size and j < coords2.size:
        # move through first list until you overlap or pass the current
        # coordinate in the second list
        while((i < coords1.size) and 
              (not ovlp(coords1[i], coords2[j])) and
              (compare(coords1[i], coords2[j]) < 0)):
            i += 1
        
        # move through second list until you overlap or pass the current
        # coordinate in the first list
        while((i < coords1.size) and (j < coords2.size) and
              (not ovlp(coords2[j], coords1[i])) and
              (compare(coords2[j], coords1[i]) < 0)):
            j += 1
        
        # keep adding overlapping coords, but keep track of where we
        # were in read list, because same coord can overlap with many
        j_overlap = j        
        while((j_overlap < coords2.size) and (i < coords1.size) and 
              ((compare(coords2[j_overlap], coords1[i]) < 0) or \
               ovlp(coords2[j_overlap], coords1[i]))):

            if ovlp(coords2[j_overlap], coords1[i]):
                overlap_list[i].append(coords2[j_overlap])

            j_overlap += 1
            
        i += 1

    # convert lists to numpy arrays
    for i in range(coords1.size):
        overlap_list[i] = np.array(overlap_list[i], dtype=coords2.dtype)

    return overlap_list



def group_by_start_end(coords):
    """Returns two dictionaries containing lists of Coords keyed by their
    start or end positions, respectively"""
    by_start = {}
    by_end = {}

    for coord in coords:
        if coord.start in by_start:
            by_start[coord.start].append(coord)
        else:
            by_start[coord.start] = [coord]
            
        if coord.end in by_end:
            by_end[coord.end].append(coord)
        else:
            by_end[coord.end] = [coord]

    return by_start, by_end



def read_bed(path, chrom_dict, min_region_size=None,
             add_one_to_start=True, other_attrib=[], 
             has_header=False):
    """Reads a list of coordinates from a BED-like file with the
    provided path (may be gzipped). If a minimum region size is
    specified then small regions are expanded symmetrically so that
    they meet this size. Unless other_attrib is specified, then only
    the first three columns of the file are used (as chromosome,
    start, end). If names of other attributes are specified, then
    these are read (as strings) and set as attributes on the returned
    coord objects"""    
    regions = []

    if path.endswith(".gz"):
        f = gzip.open(path, "rb")
    else:
        f = open(path, "r")
    
    f = open(path, "r")

    if has_header:
        header = f.readline()

    for l in f:
        line = l.rstrip()
        if line.startswith("#") or len(line) == 0:
            # skip comment and empty lines
            continue
        
        words = line.rstrip().split()

        if len(words) < 3:
            raise CoordError("BED line does not contain at "
                             "least 3 tokens:\n'%s'"
                             % line)
        
        chrom_name = words[0]
        chrom = chrom_dict[chrom_name]
        if add_one_to_start:
            start = int(words[1]) + 1
        else:
            start = int(words[1])
        end = int(words[2])

        if start < 1:
            raise CoordError("start of region (%d) is before "
                             "start of chromsome\n" % start)
        if end > chrom.length:
            raise CoordError("end of region (%d) is greater than "
                             " chromosome end (%d)" % (end, chrom.length))
        
        if min_region_size and ((end - start + 1) < min_region_size):
            # expand this region, it is less than the minimum size
            midpoint = (start + end) / 2
            half_min_size = min_region_size / 2

            start = midpoint - half_min_size
            end = midpoint + half_min_size
            
            if start < 1:
                start = 1
            if end > chrom.length:
                end = chrom.length

        region = Coord(chrom, start, end, strand=0)

        # set additional attributes on the coord object by reading 
        # them from subsequent columns on the line
        idx = 3
        for attrib_name in other_attrib:
            if idx >= len(words):
                raise CoordError("attribute '%s' could not be read "
                                 "from bed-like file, because line "
                                 "has only %d columns" % 
                                 (attrib_name, len(words)))

            if hasattr(region, attrib_name):
                raise CoordError("cannot add attribute %s to coordinate "
                                 "because this attribute already exists" %
                                 attrib_name)

            setattr(region, attrib_name, words[idx])

            idx += 1                                 
        
        regions.append(region)

    f.close()

    return regions




if __name__ == "__main__":
    chrom = Chromosome(1, "chr1", "123456789")
    
    coord1 = Coord(chrom, 1, 100, 1)
    coord2 = Coord(chrom, 80, 120, -1)
    coord3 = Coord(chrom, 90, 120, 1)
    coord4 = Coord(chrom, 400, 500, -1)

    coords = [coord1, coord2, coord3, coord4]

    print("NO-STRAND SORTING:")
    sort_coords(coords, use_strand=False)
    for c in coords:
        print("  " + str(c))

    print("overlaps:")
    overlaps = get_overlaps(coords, coords)
    for coord, ov_list in zip(coords, overlaps):
        print("  " + str(coord) + ":")
        for ov in ov_list:
            print("    " + str(ov))

    print("STRAND SORTING:")
    sort_coords(coords, use_strand=True)
    for c in coords:
        print("  " + str(c))

    print("overlaps:")
    overlaps = get_overlaps(coords, coords, use_strand=True)
    for coord, ov_list in zip(coords, overlaps):
        print("  " + str(coord) + ":")
        for ov in ov_list:
            print("    " + str(ov))


    # test with a numpy array
    coord_dtype = np.dtype([('chromosome_id', np.int16),
                            ('start', np.int32),
                            ('end', np.int32),
                            ('strand', np.int8)])

    coords1 = np.array([(1, 1, 100, 1),
                        (1, 80, 120, -1),
                        (1, 90, 120, 1),
                        (1, 400, 500, -1)], dtype=coord_dtype)

    coords2 = coords1.copy()

    print("\nNUMPY ARRAY COORDS:")
    print("NO-STRAND SORTING")
    coords1.sort(order=('chromosome_id', 'start'))
    coords2.sort(order=('chromosome_id', 'start'))
    overlaps = get_np_overlaps(coords1, coords2, use_strand=False)
    for i in range(coords1.size):
        print("%s:\n  %s"  % (str(coords1[i]), str(overlaps[i])))
    
    print("\nNUMPY ARRAY COORDS:")
    print("STRAND SORTING")
    coords1.sort(order=('chromosome_id', 'strand', 'start'))
    coords2.sort(order=('chromosome_id', 'strand', 'start'))
    overlaps = get_np_overlaps(coords1, coords2, use_strand=True)
    for i in range(coords1.size):
        print("%s:\n  %s"  % (str(coords1[i]), str(overlaps[i])))


    

