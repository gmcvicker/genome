import sys

import pysam
import tables
import argparse
import numpy as np

import genome.db


MIN_MAP_QUAL = 10
MAX_VAL = 65535

# CIGAR codes
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8


def create_carray(track, chrom):
    atom = tables.UInt16Atom(dflt=0)

    zlib_filter = tables.Filters(complevel=1, complib="zlib")

    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray


def add_read_depths(chrom, fwd_array, rev_array, bam_filename):
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

        # read pos refers to position of first mapping base, need to skip
        # cigar ops not part of alignment or inserted wrt the reference
        begin_match = False
        offset = 0
        for op, oplen in read.cigar:
            if op in [BAM_CINS, BAM_CPAD, BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
                continue
            elif op in [BAM_CMATCH, BAM_CDEL, BAM_CEQUAL, BAM_CDIFF]:
                begin_match = True
                match_start = start - 1 + offset
                match_end = start - 1 + offset + oplen
                if read.is_reverse:
                    rev_array[match_start:match_end] += 1
                else:
                    fwd_array[match_start:match_end] += 1
            if begin_match:
                offset += oplen


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default="hg18")

    parser.add_argument("--rev_track", metavar="REV_TRACK", action="store",
                        default=None,
                        help="if specified, reverse fragment read depths "
                        "are stored in REV_TRACK and forward midpoints "
                        "are stored in TRACK; otherwise all midpoints are "
                        "stored in TRACK")

    parser.add_argument("track", action="store",
                        metavar="TRACK",
                        help="name of track to store read depths in")

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
    fwd_track = gdb.create_track(args.track)

    if args.rev_track:
        rev_track = gdb.create_track(args.rev_track)
    else:
        rev_track = None

    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)

        fwd_carray = create_carray(fwd_track, chrom)

        if rev_track:
            rev_carray = create_carray(rev_track, chrom)

        # use uint32 even though we will store values as uint16s
        # to allow for possible overflows
        fwd_array = np.zeros(chrom.length, np.uint32)

        if rev_track:
            rev_array = np.zeros(chrom.length, np.uint32)
        else:
            rev_array = fwd_array

        # fill with values
        for bam_filename in args.bam_filename:
            sys.stderr.write("  %s\n  " % bam_filename)

            add_read_depths(chrom, fwd_array, rev_array, bam_filename)

            sys.stderr.write("\n")

        # threshold values to avoid integer overflow when we store them
        threshold_large_vals(fwd_array)
        if rev_track:
            threshold_large_vals(rev_array)

        fwd_carray[:] = fwd_array

        if rev_track:
            rev_carray[:] = rev_array

    fwd_track.close()

    if rev_track:
        rev_track.close()


main()
