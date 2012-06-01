import sys
import os

import pysam
import tables
import argparse
import numpy as np

import genome.db


USE_WEIGHTS = False

NUCSOME_SIZE = 147
NUCSOME_MID = 74

# from paired end reads mean frag len is 150, so mean midpoint is 75
MEAN_FRAG_MID = 75
MIN_MAP_QUAL = 10
MAX_VAL = 255

def create_carray(track, chrom):
    if USE_WEIGHTS:
        atom = tables.Float32Atom(dflt=0.0)
    else:
        atom = tables.UInt8Atom(dflt=0)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray



def add_pe_mnase_mids(chrom, fwd_array, rev_array, bam_filename,
                      min_frag_len, max_frag_len,
                      fwd_dup_counts, rev_dup_counts, max_dups):
    samfile = pysam.Samfile(bam_filename, "rb")

    count = 0
    
    for read in samfile.fetch(reference=chrom.name,
                              start=1, end=chrom.length):

        count += 1
        if count > 100000:
            sys.stderr.write(".")
            count = 0
        
        # reads appear twice, once for each side, only want to consider once
        if not (read.is_read1 or read.is_read2):
            continue
        if read.is_read2:
            continue

        # require that both sides are uniquely mapped
        if read.is_unmapped or read.mate_is_unmapped:
            continue

        if read.is_reverse == read.mate_is_reverse:
            # reads mapped to same strand...
            continue

        if read.mapq < MIN_MAP_QUAL:
            # read has poor mapping quality
            continue

        # remember pysam pos starts at 0, not 1
        if read.is_reverse:
            isize = -read.isize
            pos = read.mpos + 1
        else:
            isize = read.isize
            pos = read.pos + 1
        
        if isize < min_frag_len or isize > max_frag_len:
            continue

        if read.is_reverse:
            dup_counts = rev_dup_counts
            array = rev_array
        else:
            dup_counts = fwd_dup_counts
            array = fwd_array

        if dup_counts:
            # check if this read is a duplicate
            key = "%d:%d" % (pos, isize)
            if key in dup_counts:
                dup_counts[key] += 1
            else:
                dup_counts[key] = 1

            if dup_counts[key] > max_dups:
                # skip this read, too many duplicates
                continue
        
        # read looks good...
        # how different is this from expected size?
        size_diff = abs(NUCSOME_SIZE - isize)

        if USE_WEIGHTS:
            # assume that real midpoint occurs with uniform probability
            # within range of possible 147 bp midpoints
            mid_start = NUCSOME_MID - size_diff
            mid_end   = NUCSOME_MID + size_diff
            weight = 1.0 / (mid_end - mid_start + 1.0)
        else:
            mid_start = isize / 2
            mid_end = mid_start
            weight = 1

        # also weight by mapping error probability?
        # weight *= 10.0**(-0.1 * read.mapq)

        # convert to chromosome array index (first base on chr is 0)
        start_idx = mid_start + pos - 1
        # end idx is 'half-open' in style of python slice and UCSC coords
        end_idx = mid_end + pos

        array[start_idx:end_idx] += weight




def add_se_mnase_mids(chrom, fwd_array, rev_array, bam_filename,
                      min_frag_len, max_frag_len, 
                      fwd_dup_counts, rev_dup_counts, max_dups):
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
        # mean fragent midpoint from PE reads is 75 bp
        if read.is_reverse:
            dyad_pos = read.pos + read.alen - MEAN_FRAG_MID + 1
        else:
            dyad_pos = read.pos + MEAN_FRAG_MID

        if dyad_pos < 1:
            dyad_pos = 1
        if dyad_pos > chrom.length:
            dyad_pos = chrom.length


        if read.is_reverse:
            dup_counts = rev_dup_counts
            array = rev_array
        else:
            dup_counts = fwd_dup_counts
            array = fwd_array
            
        if dup_counts:
            # check if this read is a duplicate
            key = "%d" % pos
            if key in dup_counts:
                dup_counts[key] += 1
            else:
                dup_counts[key] = 1

            if dup_counts[key] > max_dups:
                # skip this read, too many duplicates
                continue

        arra[dyad_pos-1] += 1
        



def parse_args():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    
    
    group.add_argument("--single_end", action="store_true",
                       help="reads are single-end")
    
    group.add_argument("--paired_end", action="store_true",
                       help="reads are paired-end")

    parser.add_argument("--chrom", default=None, metavar="CHROMOSOME",
                        help="only import reads for specified chromosome or "
                        "range of chromosomes")

    parser.add_argument('--assembly', help="assembly that reads "
                       "were mapped to (e.g. hg18)",
                       default="hg18")
    
    parser.add_argument("--min_frag_size", action="store",
                        type=int, default=117,
                        help="minimum size of fragments to store")

    parser.add_argument("--max_frag_size", action="store",
                        type=int, default=172,
                        help="maximum size of fragments to store")

    parser.add_argument("--rev_track", metavar="REV_TRACK", action="store",
                        default=None,
                        help="if specified, reverse fragment midpoints "
                        "are stored in REV_TRACK and forward midpoints "
                        "are stored in TRACK; otherwise all midpoints are "
                        "stored in TRACK")

    parser.add_argument("--max_duplicates", action="store", type=int,
                        default=None, help="maximum number of duplicate "
                        "fragments (same strand, start, fragment "
                        "size) to store midpoints for")
    
    parser.add_argument("track", action="store",
                        metavar="TRACK",
                        help="name of track to store midpoints in")
    
    parser.add_argument("bam_filename", action="store", nargs="+",
                        help="sorted BAM file to read data from")


    args = parser.parse_args()

    if args.min_frag_size < 1:
        raise ValueError("--min_frag_size argument must be >= 1")

    if args.min_frag_size > args.max_frag_size:
        raise ValueError("min_frag_size must be > max_frag_size")

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

    fwd_track = gdb.create_track(args.track)

    if args.rev_track:
        rev_track = gdb.create_track(args.rev_track)
    else:
        rev_track = None

    if args.chrom:
        chromosomes = gdb.get_chromosomes_from_args([args.chrom])
    else:
        chromosomes = gdb.get_chromosomes()
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)

        if args.max_duplicates:
            fwd_dup_counts = {}
            rev_dup_counts = {}
        else:
            fwd_dup_counts = None
            rev_dup_counts = None
        
        fwd_carray = create_carray(fwd_track, chrom)

        if rev_track:
            rev_carray = create_carray(rev_track, chrom)

        if USE_WEIGHTS:
            fwd_array = np.zeros(chrom.length, np.float32)
            if rev_track:
                rev_array = np.zeros(chrom.length, np.float32)
            else:
                rev_array = fwd_array
        else:
            # use uint16 instead of uint8 to allow for possible overflows
            fwd_array = np.zeros(chrom.length, np.uint16)

            if rev_track:
                rev_array = np.zeros(chrom.length, np.uint16)
            else:
                rev_array = fwd_array
        
        # fill with values
        for bam_filename in args.bam_filename:
            sys.stderr.write("  %s\n  " % bam_filename)

            if args.paired_end:
                add_pe_mnase_mids(chrom, fwd_array, rev_array, bam_filename,
                                  args.min_frag_size, args.max_frag_size,
                                  fwd_dup_counts, rev_dup_counts,
                                  args.max_duplicates)
            else:
                add_se_mnase_mids(chrom, fwd_array, rev_array, bam_filename,
                                  args.min_frag_size, args.max_frag_size,
                                  fwd_dup_counts, rev_dup_counts,
                                  args.max_duplicates)
                
            sys.stderr.write("\n")

        if not USE_WEIGHTS:
            threshold_large_vals(fwd_array)
            if rev_track:
                threshold_large_vals(rev_array)
        
        fwd_carray[:] = fwd_array

        if rev_track:
            rev_carray[:] = rev_array
            sys.stderr.write("  stored %d fwd midpoints\n" %
                             np.sum(fwd_array))
            sys.stderr.write("  stored %d rev midpoints\n" %
                             np.sum(rev_array))
        else:
            sys.stderr.write("  stored %d midpoints\n" % np.sum(fwd_array))
    
    fwd_track.close()

    if rev_track:
        rev_track.close()


main()
        
        
    
