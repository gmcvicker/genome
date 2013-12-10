import sys
import os

import gzip
import tables
import argparse
import numpy as np

import genome.coord
import genome.db


MAX_VAL = 255



def main():
    args = parse_args()
    
    # create a database track
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    track = gdb.create_track(args.track)
    
    chromosomes = gdb.get_all_chromosomes()
    chrom_dict = gdb.get_chromosome_dict()
    
    carrays = {}
    arrays = {}

    sys.stderr.write("using frag size range: %d-%d\n"
                     % (args.min_frag_size, args.max_frag_size))

    # initialize arrays for all chromosomes
    sys.stderr.write("initializing track\n")
    for chrom in chromosomes:
        carrays[chrom.name] = create_carray(track, chrom)
        arrays[chrom.name] = np.zeros(chrom.length, dtype=np.uint8)

    # parse midpoints from mapped SOLiD reads
    for filename in args.solid_filename:
        sys.stderr.write("parsing SOLiD file %s\n" % filename)
        count_midpoints(filename, arrays, chrom_dict, 
                        args.min_frag_size, args.max_frag_size)

    sys.stderr.write("storing data in track\n")
    for chrom in chromosomes:
        carrays[chrom.name][:] = arrays[chrom.name]

    sys.stderr.write("done")
    track.close()




def create_carray(track, chrom):
    atom = tables.UInt8Atom(dflt=0)
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray




def count_midpoints(filename, chrom_arrays, chrom_dict,
                    min_frag_size, max_frag_size):
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    # skip header line
    f.readline()

    is_first = True
    read1 = read2 = None

    count = 0
    total_count = 0

    for l in f:        
        words = l.rstrip().split()

        # this is either a string that gives number of matches 
        # or a string that gives Y / N (for unique / not-unique)
        n_match_str = words[8]
        if n_match_str == 'N':
            continue
        elif n_match_str == 'Y':
            pass
        else:
            n_matches = int(n_match_str)
            if n_matches > 1:
                # skip read pairs that mapped to multiple locations
                continue

        chrom_name = "chr" + words[2]

        if chrom_name in chrom_dict:
            chrom = chrom_dict[chrom_name]
        else:
            if "mito" in chrom_name:
                chrom = chrom_dict['chrM']
            else:
                sys.stderr.write("WARNING: unknown chromosome %s\n" % 
                                 chrom.name)
                continue

        strand = genome.coord.parse_strand(words[6])
        coord = genome.coord.Coord(chrom, int(words[3]), int(words[4]),
                                   strand=strand)

        if is_first:
            # this is first half of read pair
            read1 = coord
            is_first = False
            continue

        # this is second half of read pair
        read2 = coord
        is_first = True

        if read1.chrom.name != read2.chrom.name:
            sys.stderr.write("WARNING: skipping read pair on "
                             "different chromosomes\n")
            continue

        if read1.strand == read2.strand:
            sys.stderr.write("WARNING: skipping read pair on "
                             "same strand\n")
            continue

        if read1.strand == 1:
            frag_size = read2.end - read1.start + 1
            midpoint = read1.start + (frag_size / 2)
        else:
            frag_size = read1.end - read2.start
            midpoint = read2.start + (frag_size / 2)

        if frag_size < 0:
            sys.stderr.write("WARNING: skipping read pair that implies "
                             "negative fragment size\n")
            continue

        if (frag_size < min_frag_size) or (frag_size > max_frag_size):
            # fragment size is outside acceptable range
            continue

        carray = chrom_arrays[read1.chrom.name]

        if midpoint < 1 or midpoint > read1.chrom.length:
            sys.stderr.write("WARNING: fragment midpoint %d is outside "
                             "of %s range 1-%d\n" % (midpoint,
                                                   read1.chrom.name,
                                                   read1.chrom.length))

        if carray[midpoint-1] < MAX_VAL:
            # count this midpoint
            count += 1
            total_count += 1
            if (count % 100000) == 0:
                sys.stderr.write(".")
                count = 0
            
            carray[midpoint-1] += 1
        else:
            pass
            # sys.stderr.write("WARNING: midpoint count at %s:%d "
            #                  "exceeds maximum value of %d\n"
            #                  % (read1.chrom.name, midpoint, MAX_VAL))
    
    f.close()
    sys.stderr.write("\n")
    sys.stderr.write("total_count: %d\n" % total_count)

    for chrom_name in chrom_arrays.keys():
        carray = chrom_arrays[chrom_name]
        sys.stderr.write("%s: sites with max val (%d): %d\n" % 
                         (chrom_name, MAX_VAL, np.sum(carray[:] == MAX_VAL)))


    
    
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--assembly', help="assembly that reads "
                       "were mapped to (e.g. hg18)",
                       default="hg18")
    
    parser.add_argument("--min_frag_size", action="store",
                        type=int, default=117,
                        help="minimum size of fragments to store")

    parser.add_argument("--max_frag_size", action="store",
                        type=int, default=172,
                        help="maximum size of fragments to store")
    
    parser.add_argument("track", action="store",
                        metavar="TRACK",
                        help="name of track to store midpoints in")
    
    parser.add_argument("solid_filename", action="store", nargs="+",
                        help="file containing mapped SOLiD reads")

    args = parser.parse_args()

    if args.min_frag_size < 1:
        raise ValueError("--min_frag_size argument must be >= 1")

    if args.min_frag_size > args.max_frag_size:
        raise ValueError("min_frag_size must be > max_frag_size")

    return args
    






main()
        
        
    
