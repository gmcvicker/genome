#!/bin/env python

"""This program lifts over data contained in an HDF5 file from one
assembly to another using a mapping defined by a UCSC chain
file. Currently this program only works for 1D arrays. Regions where
no data were mapped to in the new assembly are set to nan for floating
point data types, or 0 for integer datatypes. This program currently
takes a single un-stranded track as input, but it could be easily
modified to take a pair of fwd/rev strand tracks instead.
"""

import sys
import argparse
import gzip

import tables

import genome.db
import genome.coord
from genome.coord import Coord, CoordError


class Chain(object):
    def __init__(self, line, from_chrom_dict, to_chrom_dict):
        words = line.rstrip().split()

        #
        # Parse from coordinates for the this chain, being careful
        # that they look correct
        #
        from_chrom_name = words[2]
        from_chrom_size = int(words[3])

        if from_chrom_name in from_chrom_dict:
            from_chrom = from_chrom_dict[from_chrom_name]
        else:
            raise CoordError("chromosome '%s' does not exist "
                             "in 'from' database'" % from_chrom_name)

        if from_chrom.length != from_chrom_size:
            raise CoordError("chromosome length mismatch between chain "
                             "file and 'from' database %s:1-%d != %s:1-%d"
                             % (from_chrom_name, from_chrom_size,
                                from_chrom_name, from_chrom.length))

        if words[4] == '+':
            from_strand = 1
        else:
            raise CoordError("expected from strand from be '+'")
        
        from_start = int(words[5]) + 1
        from_end   = int(words[6])
        self.from_coord = Coord(from_chrom, from_start, from_end,
                                strand=from_strand)

        #
        # now parse "to" coordaintes
        #
        to_chrom_name = words[7]                
        to_chrom_size = int(words[8])

        if to_chrom_name in to_chrom_dict:
            to_chrom = to_chrom_dict[to_chrom_name]
        else:
            raise CoordError("chromosome %s does not exist "
                             "in 'to' database'" % to_chrom_name)

        if to_chrom.length != to_chrom_size:
            raise CoordError("chromosome length mismatch between chain "
                             "file and 'to' database %s:1-%d != %s1-%d"
                             % (to_chrom_name, to_chrom_size,
                                to_chrom_name, to_chrom.length))

        if words[9] == '+':
            to_strand = 1
        elif words[9] == '-':
            to_strand = -1
        else:
            raise CoordError("expected 'to' strand to be '+' or '-', "
                             "not '%s'" % words[9])

        if to_strand == 1:
            to_start = int(words[10]) + 1
            to_end = int(words[11])
        else:
            # coordinates in file are on reverse strand
            # and are 0-based half-open
            to_start = to_chrom_size - int(words[11]) + 1
            to_end = to_chrom_size - int(words[10])

        self.to_coord = Coord(to_chrom, to_start,
                              to_end, strand=to_strand)

        self.ori = to_strand

        
        
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("from_assembly", metavar="FROM_ASSEMBLY",
                        help="Assembly to copy track FROM (e.g. hg18)")

    parser.add_argument("to_assembly", metavar="TO_ASSEMBLY",
                        help="Assembly to copy track TO (e.g. hg19)")

    parser.add_argument("liftover_file", metavar="LIFTOVER_FILE",
                        help="Path to UCSC chain file that "
                        "describes the coordinate conversion between "
                        "assemblies (e.g. hg18ToHg19.over.chain.gz)")
    
    parser.add_argument("track", metavar="TRACK", help="Name of track we want "
                        "to convert from old assembly to new assembly")

    parser.add_argument("--rev_track", metavar="REV_TRACK",
                        default=None,
                        help="If this is argument is specified then two "
                        "tracks are lifted over at the same time--one for each "
                        "strand. In regions where the assemblies differ in "
                        "orientation, the data from TRACK and REV_TRACK are "
                        "swapped")

    args = parser.parse_args()

    return args




def create_carray(from_track, from_chrom, to_track, to_chrom):
    """Creates a new CArray node on the 'to' track, copying the
    datatype of the provided 'from' track."""

    from_node = from_track.get_array(from_chrom)

    if len(from_node.shape) != 1:
        raise NotImplementedError("only liftover of 1D arrays is implemented")

    if from_node.shape[0] != from_chrom.length:
        raise CoordError("expected length of array (%d) to match "
                         "chromosome length (%d)" %
                         (from_node.shape[0], from_chrom.length))

    atom = from_node.atom.copy()

    zlib_filter = tables.Filters(complevel=1, complib="zlib")

    if str(atom.dtype).startswith("int") or str(atom.dtype).startswith('uint'):
        atom.dflt = 0
    elif str(atom.dtype).startswith("float"):
        atom.dflt = np.nan
    else:
        raise ValueError("unknown datatype '%s', expected int or float" %
                         str(atom.dtype))
    
    # create new CArray for this chromosome
    shape = [to_chrom.length]
    carray = to_track.h5f.createCArray(to_track.h5f.root, to_chrom.name,
                                       atom, shape, filters=zlib_filter)

    
    return carray




def copy_data(from_track, rev_from_track, from_chrom, from_start, from_end,
              to_track, rev_to_track, to_chrom, to_start, to_end, ori):
    """Copy data from one track to another using the provided
    coordinates."""

    sys.stderr.write("%s:%d-%d => %s:%d-%d(%d)\n" %
                     (from_chrom.name, from_start, from_end,
                      to_chrom.name, to_start, to_end, ori))

    if not from_track.has_chromosome(from_chrom):
        # ignore this region, was not present in original track
        return

    if to_track.has_chromosome(to_chrom):
        to_array = to_track.get_array(to_chrom)

        if rev_to_track:
            rev_to_array = rev_to_track.get_array(to_chrom)
    else:
        # need to create this region on new tracks
        to_array = create_carray(from_track, from_chrom, to_track, to_chrom)

        if rev_to_track:
            rev_to_array = create_carray(rev_from_track, from_chrom,
                                         rev_to_track, to_chrom)
    
    # retrieve values from 'from' track
    from_vals = from_track.get_nparray(from_chrom, from_start, from_end)

    if rev_to_track:
        rev_from_vals = rev_from_track.get_nparray(from_chrom, from_start,
                                                   from_end)
        if ori == -1:
            # swap values from fwd/rev tracks (and reverse direction)
            tmp = rev_from_vals[::-1]
            rev_from_vals = from_vals[::-1]
            from_vals = tmp

        to_array[(to_start-1):to_end] = from_vals
        rev_to_array[(to_start-1):to_end] = rev_from_vals
    else:
        # need to flip the values if the orientation is flipped
        if ori == -1:
            from_vals = from_vals[::-1]
        
        # copy values to 'to' track
        to_array[(to_start-1):to_end] = from_vals


    


def liftover_data(chain_file, from_gdb, to_gdb, from_track, rev_from_track,
                  to_track, rev_to_track):
    # read chromosomes from both databases
    from_chrom_dict = from_gdb.get_chromosome_dict()
    to_chrom_dict = to_gdb.get_chromosome_dict()

    if chain_file.endswith(".gz"):
        f = gzip.open(chain_file)
    else:
        f = open(chain_file)

    for line in f:
        if line.startswith('chain'):
            # this is the start of a new chain
            chain = Chain(line, from_chrom_dict, to_chrom_dict)
            # set starts (or ends) of block to the edge of this chain
            from_start = chain.from_coord.start
            from_end = None
            if chain.ori == 1:
                to_start = chain.to_coord.start
                to_end = None
            else:
                to_end = chain.to_coord.end
                to_start = None
        else:
            # this is a line describing a block within a chain
            # update ends of blocks
            words = line.rstrip().split()

            if len(words) == 0:
                # there is a blank line between chain lines
                continue
            
            size = int(words[0])

            from_end = from_start + size - 1

            if chain.ori == 1:
                to_end = to_start + size - 1
            else:
                # on reverse strand: to coords go other direction
                to_start = to_end - size + 1

            # now liftover the data data for this block
            copy_data(from_track, rev_from_track, chain.from_coord.chrom,
                      from_start, from_end, to_track, rev_to_track,
                      chain.to_coord.chrom, to_start, to_end, chain.ori)
            
            if len(words) == 3:
                # there are gaps between blocks.
                # update the start of the next block to be past the end
                # of this one by a defined offset
                from_offset = int(words[1])
                to_offset = int(words[2])

                from_start = from_end + from_offset + 1

                if chain.ori == 1:
                    to_start = to_end + to_offset + 1
                else:
                    to_end = to_start - to_offset - 1
                    
            elif len(words) == 1:
                # this was the last line in chain... don't need to do anything
                # sanity check: we should be at the end of the chain
                if from_end != chain.from_coord.end:
                    raise CoordError("end of last block (%d) does not "
                                     "match end of chain: %d" %
                                     (from_end, chain.from_coord.end))
                
                if chain.ori == 1:
                    if to_end != chain.to_coord.end:
                        raise CoordError("end of last block (%d) does not "
                                         "match end of chain: %d" %
                                         (to_end, chain.to_coord.end))
                else:
                    if to_start != chain.to_coord.start:
                        raise CoordError("start of last block (%d) does not "
                                         "match start of chain: %d" %
                                         (to_start, chain.to_coord.start))
                        
            else:
                raise ValueError("expected line to have 1 or 3 tokens")

    f.close()


def main():
    args = parse_args()

    from_gdb = genome.db.GenomeDB(assembly=args.from_assembly)
    to_gdb = genome.db.GenomeDB(assembly=args.to_assembly)

    # get the original track
    from_track = from_gdb.open_track(args.track)

    # create new track with same name in new database
    to_track = to_gdb.create_track(args.track)

    if args.rev_track:
        # there are separate forward and reverse strand tracks
        rev_from_track = from_gdb.open_track(args.rev_track)
        rev_to_track = to_gdb.create_track(args.rev_track)
    else:
        rev_from_track = None
        rev_to_track = None
        
                        
    liftover_data(args.liftover_file, from_gdb, to_gdb, 
                  from_track, rev_from_track, to_track, rev_to_track)
        
    from_track.close()
    to_track.h5f.flush()
    to_track.close()

    if args.rev_track:
        rev_from_track.close()
        rev_to_track.h5f.flush()
        rev_to_track.close()
    
    sys.stderr.write("done\n")
    

main()
