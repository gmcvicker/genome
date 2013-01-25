
import argparse
import numpy as np
import sys

import scipy.stats
import genome.db
import tables

MIN_SHIFT = 20
MAX_SHIFT = 250
MAX_SHIFT_MAX_DIST = 50
SMOOTH_WIN_SIZE = 50



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default="hg18")

    parser.add_argument("--dtype", choices=("uint8", "uint16"),
                        default="uint8", help="datatype of combined track")
    
    parser.add_argument("fwd_track", action="store",
                        metavar="FWD_TRACK",
                        help="name of track to store forward read counts in")

    parser.add_argument("rev_track", metavar="REV_TRACK", action="store",
                        help="name of track to store reverse read counts in")
    
    parser.add_argument("combined_track", metavar="COMBINED_TRACK",
                        action="store", help="name of new track to store "
                        "combined counts in")
    
    args = parser.parse_args()

    return args




def create_carray(track, chrom, dtype):
    if dtype == 'uint8':
        atom = tables.UInt8Atom(dflt=0)
    elif dtype == 'uint16':
        atom = tables.UInt16Atom(dflt=0)
    else:
        raise NotImplementedError("support for dtype %s not "
                                  "yet implemented" % dtype)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray


def find_max_cov(fwd_vals, rev_vals):    
    max_cov = 0.0
    max_cov_shift = 0

    for i in range(MIN_SHIFT, MAX_SHIFT):
        # shift fwd values ahead by i bp
        # shift rev values back by i bp
        matrix = np.vstack((np.roll(fwd_vals, i),
                            np.roll(rev_vals, -i)))

        # compute covariance between fwd and rev
        cov = np.cov(matrix)[0,1]

        max_str = ""
        if cov > max_cov:
            # this is the highest covariance we've seen yet
            max_cov = cov
            max_cov_shift = i
            max_str = " *"

        last_max = i - max_cov_shift

        if last_max > MAX_SHIFT_MAX_DIST:
            # we seem to be well past the peak covariance
            break

        sys.stderr.write("shift: %d, cov: %g%s\n" % (i, cov, max_str))
        
    sys.stderr.write("shift: %d, cov: %g\n" % (max_cov_shift, max_cov))

    return max_cov_shift


def find_strand_shift(gdb, fwd_track, rev_track):
    # use a single chromosome to find shift between strands
    # that gives max covariance.
    if gdb.assembly == "dm3":
        chrom = gdb.get_chromosome("chr2L")
    else:
        chrom = gdb.get_chromosome("chr21")

    sys.stderr.write("retrieving values\n")
    fwd_vals = fwd_track.get_nparray(chrom)
    rev_vals = rev_track.get_nparray(chrom)

    # smooth values using sliding window
    sys.stderr.write("smoothing values\n")
    win = np.ones(SMOOTH_WIN_SIZE)
    fwd_vals = np.convolve(fwd_vals, win, mode="same")
    rev_vals = np.convolve(rev_vals, win, mode="same")

    # only use regions with high density of sites so all of contribution
    # to covariance comes from these
    fwd_cutoff = scipy.stats.scoreatpercentile(fwd_vals, 95)
    rev_cutoff = scipy.stats.scoreatpercentile(rev_vals, 95)
    fwd_vals[fwd_vals < fwd_cutoff] = 0
    rev_vals[rev_vals < rev_cutoff] = 0

    # find the shift that yields the max covariance
    max_cov_shift = find_max_cov(fwd_vals, rev_vals)

    return max_cov_shift
    


def main():
    args = parse_args()

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    fwd_track = gdb.open_track(args.fwd_track)
    rev_track = gdb.open_track(args.rev_track)
    combined_track = gdb.create_track(args.combined_track)
    
    shift = find_strand_shift(gdb, fwd_track, rev_track)

    sys.stderr.write("shifting fwd/rev by +%d/-%d bp\n" % (shift, shift))

    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)
        carray = create_carray(combined_track, chrom, args.dtype)

        fwd_vals = fwd_track.get_nparray(chrom)
        rev_vals = rev_track.get_nparray(chrom)

        # shift fwd / rev values by the offset that gave
        # the maximum covariance2
        fwd_vals = np.roll(fwd_vals, shift)
        rev_vals = np.roll(rev_vals, -shift)
        fwd_vals[:shift] = 0
        rev_vals[-shift:] = 0

        carray[:] = fwd_vals + rev_vals
        
    fwd_track.close()
    rev_track.close()
    combined_track.close()

    


main()
