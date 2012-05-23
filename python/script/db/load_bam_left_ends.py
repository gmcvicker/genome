
import sys
import os

import pysam
import tables
import argparse
import numpy as np

import genome.db


MIN_MAP_QUAL = 10
MAX_VAL = 65535

def create_carray(track, chrom):
    atom = tables.UInt16Atom(dflt=0)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray



def add_read_counts(chrom, array, bam_filename):
    samfile = pysam.Samfile(bam_filename, "rb")

    count = 0
    
    for read in samfile.fetch(reference=chrom.name, start=1,
                              end=chrom.length):
        count += 1
        if count > 100000:
            sys.stderr.write(".")
            count = 0

        if read.is_unmapped:
            continue

        if read.mapq < MIN_MAP_QUAL:
            # read has poor mapping quality
            continue

        # remember pysam pos starts at 0, not 1
        start = read.pos + 1
        end = start + read.alen - 1

        if end > chrom.length:
            sys.stderr.write("  warning: read coordinates %d-%d fall outside "
                             "%s range 1-%d" % (start, end, chrom.name,
                                                 chrom.length))
            end = chrom.length
        
        if start < 1:
            sys.stderr.write("  warning: read coordinates %d-%d fall outside "
                             "%s range 1-%d" % (start, end, chrom.name,
                                                 chrom.length))
            start = 1

        array[start-1] += 1



def parse_args():    
    parser = argparse.ArgumentParser(description="Stores the number of "
                                     "reads that map to each position "
                                     "in the genome. Ths script ignores "
                                     "strand and defines the genomic "
                                     "position of the read as the left "
                                     "end of where it mapped")

    parser.add_argument("--assembly", metavar="ASSEMBLY",
                        help="genome assembly that reads were mapped "
                        "to (e.g. hg18)", default='hg18')
    
    parser.add_argument("track", action="store", metavar="TRACK",
                        help="name of track to store read counts in")
    
    parser.add_argument("bam_filename", action="store", nargs="+",
                        help="sorted BAM file to read data from")

    args = parser.parse_args()

    return args
    



def threshold_large_vals(array):
    n_large_vals = np.sum(array > MAX_VAL)
    if n_large_vals > 0:
        sys.stderr.write("%d sites exceed max value %d\n" %
                         (n_large_vals, MAX_VAL))
        array[array > MAX_VAL] = MAX_VAL



def main():
    args = parse_args()
    
    # create a database track
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    track = gdb.create_track(args.track)
    
    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)
        
        carray = create_carray(track, chrom)

        # use uint32 even though we will store values as uint16s
        # to allow for possible overflows
        array = np.zeros(chrom.length, np.uint32)
        
        # fill with values
        for bam_filename in args.bam_filename:
            sys.stderr.write("  %s\n  " % bam_filename)
            add_read_counts(chrom, array, bam_filename)
            sys.stderr.write("\n")

        # threshold values to avoid integer overflow when we store them
        threshold_large_vals(array)
        
        carray[:] = array
    
    track.close()


main()
        
        
    
