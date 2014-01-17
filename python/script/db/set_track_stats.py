#!/usr/bin/python

import sys
import numpy as np
import argparse

import genome.db
import genome.trackstat as trackstat


def parse_args():
    parser = argparse.ArgumentParser(description="Saves a number of "
                                     "statistics (max, min, mean, etc.) "
                                     "as track attributes. These can be "
                                     "rapidly retrieved by other programs.")

    parser.add_argument("--assembly", metavar="ASSEMBLY", default=None,
                        help="name of assembly (e.g. hg18)")

    parser.add_argument("track_name", metavar="TRACK_NAME",
                        help="name of track")
    
    args = parser.parse_args()

    return args


def main():    
    args = parse_args()

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    track = gdb.open_track(args.track_name, "a")
    track_stat = trackstat.set_stats(gdb, track)
    sys.stderr.write("combined %s\n" % str(track_stat))
    track.close()


if __name__ == "__main__":
    main()
