#!/bin/env python

"""Makes 2D arrays in HDF5 format that can be used to quickly convert
coordinates from one assembly (the 'from' assembly)' to another (the
'to' assembly). Each chromosome in the 'from' assembly is represented
by an N x 3 array where N is the number of bases in the chromosome.
Each row corresponds to a base position on the 'from' chromosome and
the three values give the equivalent coordinates in the 'to'
assembly. The three values are: chromosome_id, chromosome_position,
orientation (1 if same, -1 if reversed).
  
[-1,-1,0] are used to flag positions on the 'from' assembly that do not map 
to the 'to' assembly.

"""

import sys
import argparse
import gzip
import numpy as np

import tables

import genome.db
import chain

POS_UNDEF = -1

ZLIB_FILTER = tables.Filters(complevel=1, complib='zlib')

        
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
    
    parser.add_argument("tabname", metavar="TABLE", help="Name of table to store "
                        "liftover information in")

    args = parser.parse_args()

    return args




def create_tables(from_gdb, tab_name):
    track = from_gdb.create_track(tab_name)

    chrom_list = from_gdb.get_all_chromosomes()

    atom = tables.Int32Atom(dflt=POS_UNDEF)
    
    for chrom in chrom_list:        
        sys.stderr.write(" %s\n" % chrom.name)
        shape = (chrom.length, 3)
        carray = track.h5f.createCArray(track.h5f.root, chrom.name, atom, shape,
                                        filters=ZLIB_FILTER)

        #carray[:,:] = POS_UNDEF

    return track
    


def main():
    args = parse_args()

    from_gdb = genome.db.GenomeDB(assembly=args.from_assembly)
    to_gdb = genome.db.GenomeDB(assembly=args.to_assembly)

    chain_list = chain.read_chain_file(args.liftover_file, from_gdb, to_gdb)

    sys.stderr.write("initializing tables\n")
    track = create_tables(from_gdb, args.tabname)

    sys.stderr.write("writing to tables\n")
    for c in chain_list:
        sys.stderr.write("  current chain: %s\n" % str(c))
        
        tab = track.h5f.getNode("/%s" % c.from_coord.chrom.name)
        to_chr_id = c.to_coord.chrom.idnum

        for block in c.blocks:
            # update table with mappings to new chromosomal locations
            to_pos = np.arange(block.to_start, block.to_end+1)
            block_size = block.to_end - block.to_start + 1
            to_chr_id = np.repeat(c.to_coord.chrom.idnum, block_size)
            ori = np.repeat(c.ori, block_size)
            vals = np.transpose(np.vstack((to_chr_id, to_pos, ori)))
            tab[block.from_start-1:block.from_end, :] = vals
    
    

main()
