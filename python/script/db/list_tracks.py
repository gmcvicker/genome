#!/usr/bin/python

import argparse
import genome.db


parser = argparse.ArgumentParser()

parser.add_argument("--assembly", metavar="ASSEMBLY", default=None,
                    help="Assembly to list tracks for (e.g. hg18)")

parser.add_argument("--paths", dest="paths", action="store_true",
                    default=False, help="print full paths to "
                    ".h5 files, rather than tracknames")

                    
args = parser.parse_args()

gdb = genome.db.GenomeDB(assembly=args.assembly)

tracknames = gdb.list_tracks()
tracknames.sort()

for trackname in tracknames:
    if args.paths:
        print(gdb.get_track_path(trackname))
    else:
        print(trackname)
