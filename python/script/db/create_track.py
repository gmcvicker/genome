#!/bin/env python
import sys, re

from argparse import ArgumentParser

import numpy as np
import tables

import genome.track
import genome.chrom

# trackreader contains Cython bindings to C library for
# speedy parsing of large text files
import trackreader


# according to benchmarks in PyTables manual, zlib compression level 1
# works about as well as higher compression levels and is faster
ZLIB_FILTER = tables.Filters(complevel=1, complib='zlib')


def extract_chrom_name(filename):
    # does filename contain a random chromosome?
    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?\_random)", filename)
    if matches:
        return matches[0]

    # is this an alt-haplotype chromosome?
    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?\_\w+\_hap\d+)", filename)
    if matches:
        return matches[0]

    # is this a 'gl' chromosome
    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?\_gl\d+\_random)", 
                         filename)
    if matches:
        return matches[0]

    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?\_gl\d+)", filename)
    if matches:
        return matches[0]

    matches = re.findall(r"(chrUn\_gl\d+)", filename)
    if matches:
        return matches[0]

    

    # is this a fly heterochromatin chromosome?
    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?Het)", filename)
    if matches:
        return matches[0]

    # is this a fly 'extra' chromosome?
    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?extra)", filename)
    if matches:
        return matches[0]

    # check for vanilla chromosome name
    matches = re.findall(r"(chr[0-9UWXYM]+[A|B|L|R]?)[\.\_]?", filename)
    if matches:
        return matches[0]

    raise ValueError("could not parse chromosome name from filename '%s'" %
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
                        required=True,
                        choices=("fasta", "wiggle", "bedgraph", "txtfile"),
                        default="wiggle", help="format of input files")

    parser.add_argument("-a", "--assembly", default=None,
                        required=True,
                        help='genome assembly that track is on (e.g. hg19)')

    parser.add_argument("-c", "--chrom", required=True,
                        default=None,
                        help="path to chromInfo file containing "
                        "names and lengths of chromosomes")
    
    parser.add_argument("-p", "--pos_idx", action="store",
                        dest="pos_idx", type=int, default=-1,
                        help="index of position column in the "
                        "tab- or space-delimited input txtfile")

    parser.add_argument("-v", "--val_idx", action="store",
                        dest="val_idx", type=int, default=-1,
                        help="index of value column in the "
                        "tab- or space-delimited input txtfile")

    parser.add_argument("-n", "--name", action="store",
                        dest="name", default=None,
                        required=False, help="name of track")

    parser.add_argument("-s", "--desc", action="store",
                        dest="desc", default="",
                        required=False, help="description of track")

    parser.add_argument("hdf5_outfile", action="store",
                        help="output file to store data in")
    
    parser.add_argument("filename", action="store", nargs="+",
                        help="input file to read data from")
    
    options = parser.parse_args(args)

    if options.format == "txtfile":
        if options.pos_idx < 0 or options.val_idx < 0:
            parser.print_help()
            parser.print_usage()
            parser.error("positive pos_idx and val_idx values must be "
                         "provided in order to parse txtfiles")

    if options.format == "fasta":
        if options.dtype != "uint8":
            # only uint8 makes sense for storing DNA sequence
            options.dtype = "uint8"
            sys.stderr.write("using uint8 datatype for sequence data\n")
    

    return options



def main(options):
    chrom_dict = genome.chrom.parse_chromosomes_dict(options.chrom)
    chrom_list = list(chrom_dict.values())

    track = genome.track.Track(options.hdf5_outfile, mode="w",
                               assembly=options.assembly,
                               chromosomes=chrom_list,
                               name=options.name, desc=options.desc)

    
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
        chrom_name = extract_chrom_name(filename)
        
        if chrom_name not in chrom_dict:
            raise ValueError("unknown chromosome '%s'" % chrom_name)
        
        chrom = chrom_dict[chrom_name]
        sys.stderr.write("%s (len: %d)\n" % (chrom.name, chrom.length))

        # create a chunked array with one dimension the length
        # of the chromosome
        shape = [chrom.length]
        carray = track.h5f.create_carray(track.h5f.root, chrom_name,
                                         atom, shape, filters=ZLIB_FILTER)

        
        # populate the array with data read from a file
        carray[:] = trackreader.read_file(path, chrom,
                                          dtype=options.dtype,
                                          format=options.format,
                                          pos_idx=options.pos_idx,
                                          val_idx=options.val_idx)

    track.close()


if __name__ == "__main__":
    options = parse_options(sys.argv[1:])
    main(options)
