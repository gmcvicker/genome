#!/usr/bin/python

import argparse
import genome.db


parser = argparse.ArgumentParser()
parser.add_argument("--assembly", metavar="ASSEMBLY", default="hg18",
                    help="Assembly to list tracks for (e.g. hg18)")
args = parser.parse_args()


gdb = genome.db.GenomeDB(assembly=args.assembly)

tracknames = gdb.list_tracks()
tracknames.sort()

for trackname in tracknames:
    print trackname
