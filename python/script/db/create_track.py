#!/bin/env python

import sys, re

from argparse import ArgumentParser

# import NumPy and PyTables
import numpy as np
import tables

# import hdf5 database
import genome.db

# trackreader contains Cython bindings to C library for
# speedy parsing of large text files
import trackreader


# according to benchmarks in PyTables manual, zlib compression level 1
# works about as well as higher compression levels and is faster
ZLIB_FILTER = tables.Filters(complevel=1, complib='zlib')


def extract_chrom_name(filename):
    # does filename contain a random chromosome?
    matches = re.findall(r"(chr[0-9XYM]+\_random)", filename)
    if matches:
        return matches[0]

    # is this an alt-haplotype chromosome?
    matches = re.findall(r"(chr[0-9XYM]+\_\w+\_hap\d+)", filename)
    if matches:
        return matches[0]

    # check for vanilla chromosome name
    matches = re.findall(r"(chr[0-9XYM]+)", filename)
    if matches:
        return matches[0]

    raise ValueError("could not parse chromosome name from filename %s" %
                     filename)

    
        

def print_examples():
    progname = sys.argv[0]
    sys.stderr.write("Examples:\n")
    sys.stderr.write("  %s --dtype=int8 --format=txtfile --pos_idx=0 \\\n"
                     "     --val_idx=1 mappability_35bp "
                     "data/*mappability.gz\n\n" % progname)
    sys.stderr.write("  %s --dtype=float32 --format=wiggle \\\n"
                     "     avg_read_depth_35bp data/chr*.wig.gz\n\n" % progname)
    sys.stderr.write("  %s --dtype=int16 --format=bedgraph \\\n"
                     "     /encode/pol2_chip_seq/GM10847 data/chr*.bedgraph.gz\n\n" 
                     % progname)
    sys.stderr.write("  %s --dtype=uint8 --format=fasta \\\n"
                     "     seq data/hg18/chr*.fa.gz\n\n" 
                     % progname)


def parse_options(args):
    parser = ArgumentParser()

    parser.add_argument("-d", "--dtype", action="store", dest="dtype",
                        choices=("int8", "uint8", "int16", "float32"),
                        default="float32", help="datatype of values to store")

    parser.add_argument("-f", "--format", action="store", dest="format",
                        choices=("fasta", "wiggle", "bedgraph", "xb", 
                                 "txtfile"),
                        default="wiggle", help="format of input files")

    parser.add_argument("-p", "--pos_idx", action="store",
                        dest="pos_idx", type=int, default=-1,
                        help="index of position column in the "
                        "tab- or space-delimited input txtfile")

    parser.add_argument("-v", "--val_idx", action="store",
                        dest="val_idx", type=int, default=-1,
                        help="index of value column in the "
                        "tab- or space-delimited input txtfile")

    parser.add_argument("-s", "--strand", action="store",
                        dest="strand", default="forward",
                        choices=("forward", "reverse"),
                        help="strand of data to import (for xb files only)")

    parser.add_argument("track_name", action="store", nargs=1,
                        help="name of track to store data in")
    
    parser.add_argument("filename", action="store", nargs="+",
                        help="input file to read data from")
    
    options = parser.parse_args(args)

    if options.format == "txtfile":
        if options.pos_idx < 0 or options.val_idx < 0:
            parser.print_help()
            parser.print_usage()
            parser.error("positive pos_idx and val_idx values must be "
                         "provided in order to parse txtfiles")

    return options



def main(options):
    gdb = genome.db.GenomeDB()

    chrom_dict = gdb.get_chromosome_dict()

    track = gdb.create_track(options.track_name[0])
    
    if options.dtype == "float32":
        atom = tables.Float32Atom()
    elif options.dtype == "int8":
        atom = tables.Int8Atom()
    elif options.dtype == "uint8":
        atom = tables.UInt8Atom()
    elif options.dtype == "int16":
        atom = tables.Int16Atom()
    else:
        raise NotImplementedError("datatype %s not implemented" % dtype)

    for path in options.filename:
        filename = path.split("/")[-1]

        if options.format in ("xb", "xbf"):
            # all of the chromosomes are in a single file...
            chrom_names = [chrom.name for chrom in gdb.get_chromosomes()]
        else:
            chrom_names = [extract_chrom_name(filename)]
            
        for chrom_name in chrom_names:
            if chrom_name not in chrom_dict:
                raise ValueError("unknown chromosome '%s'" % chrom_name)

            chrom = chrom_dict[chrom_name]
            sys.stderr.write(chrom_name + "\n")

            # create a chunked array with one dimension the length
            # of the chromosome
            shape = [chrom.length]
            carray = track.h5f.createCArray(track.h5f.root, chrom_name,
                                            atom, shape, filters=ZLIB_FILTER)

            # populate the array with data read from a file
            carray[:] = trackreader.read_file(path, chrom,
                                              dtype=options.dtype,
                                              format=options.format,
                                              pos_idx=options.pos_idx,
                                              val_idx=options.val_idx,
                                              strand=options.strand)

    track.close()


if __name__ == "__main__":
    options = parse_options(sys.argv[1:])
    main(options)
