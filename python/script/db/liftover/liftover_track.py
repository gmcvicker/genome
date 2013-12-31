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

from chain import Chain, ChainBlock, read_chain_file


        
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

    # sys.stderr.write("%s:%d-%d => %s:%d-%d(%d)\n" %
    #                  (from_chrom.name, from_start, from_end,
    #                   to_chrom.name, to_start, to_end, ori))

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

    sys.stderr.write("reading chain file\n")
    chain_list = read_chain_file(chain_file, from_gdb, to_gdb)

    sys.stderr.write("copying data\n")
    for chain in chain_list:
        sys.stderr.write("  current chain: %s\n" % str(chain))
        for block in chain.blocks:
            # liftover the data data for this block
            copy_data(from_track, rev_from_track, chain.from_coord.chrom,
                      block.from_start, block.from_end, to_track, rev_to_track,
                      chain.to_coord.chrom, block.to_start, block.to_end, chain.ori)
            

            
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
    
    

main()
