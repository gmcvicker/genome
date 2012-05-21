import sys
import tables
import argparse

import genome.db


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", default='hg18',
                        help="assembly to print chromosomes for"
                        " (e.g. hg18)")

    return parser.parse_args()


def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    for chrom in gdb.get_all_chromosomes():
        print "%s\t%d" % (chrom.name, chrom.length)

main()
