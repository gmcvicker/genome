import sys
import os
import gzip
import argparse
import subprocess

import numpy as np

import genome.db
import genome.wig


    
def parse_args():
    parser = argparse.ArgumentParser(description="dumps a track to "
                                     "a wiggle file(s)")

    parser.add_argument("--assembly", action="store",
                        default="hg19", help="genome assembly to use")
    
    parser.add_argument("--chrom", action="store",
                        default=None, 
                        help="range of chromosomes to run on")

    parser.add_argument("--combine_files", action="store_const",
                        const=True, default=False,
                        help="combine chromosome files into one file "
                        "(default=no)")
    
    parser.add_argument("track_name", action="store",
                        help="name of track to read data from")

    
    parser.add_argument("output_dir", action="store",
                        default=None,
                        help="write separate files for each chromosome to "
                        "this directory, otherwise write to stdout")

    return parser.parse_args()



def combine_files(output_dir, filenames):
    filename = "%s/combined.wig.gz" % output_dir
    sys.stderr.write("combining files into file %s\n" % filename)
    out_file = open(filename, "w")

    cmd = ["gunzip", "-c"] + filenames
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=out_file)
    p1.stdout.close()
    out_file.close()
    
    

def main():
    args = parse_args()

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    out_filenames = []

    if args.chrom is None:
        # use full set of chromosomes
        chromosomes = gdb.get_chromosomes()
    else:
        # use specified chromosomes
        chromosomes = gdb.get_chromosomes_from_args(args.chrom)
        
    track = gdb.open_track(args.track_name, "r")

    
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)
        
        sys.stderr.write("  retrieving values\n")
        vals = track.get_nparray(chrom)

        # write to chromosome wiggle files
        out_filename = args.output_dir + "/%s.wig.gz" % chrom.name
        
        if os.path.exists(out_filename):
            raise IOError("output file %s already exists" % out_filename)

        out_filenames.append(out_filename)
        
        if vals.dtype == 'uint8':
            genome.wig.write_uint8(out_filename, vals, chrom.name)
        elif vals.dtype == 'float32':
            genome.wig.write_float32(out_filename, vals, chrom.name)
        else:
            raise NotImplementedError("only uint8 and float32 datatypes "
                                      "are currently implemented")

    if args.combine_files:
        combine_files(args.output_dir, out_filenames)

    track.close()



main()
