#!/usr/bin/python

import sys
import numpy as np
import argparse

import genome.db
import genome.trackstat as trackstat


def parse_args():
    parser = argparse.ArgumentParser(description="Write track "
                                     "statistics (max, min, mean, etc.) "
                                     "to stdout. Stats must already have "
                                     "been stored as attributes using "
                                     "set_track_stats.py")

    parser.add_argument("--assembly", metavar="ASSEMBLY", default=None,
                        help="name of assembly (e.g. hg18)")

    parser.add_argument("--verbose", action="store_true",
                        help="print extra info to stderr "
                        "(stats for each chromosome)")

    parser.add_argument("track_name", metavar="TRACK_NAME",
                        help="name of track")
    
    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    track = gdb.open_track(args.track_name)
    
    track_stat = trackstat.get_stats(gdb, track, verbose=args.verbose)
    sys.stdout.write("combined %s\n" % str(track_stat))

    track.close()


if __name__ == "__main__":
    main()
