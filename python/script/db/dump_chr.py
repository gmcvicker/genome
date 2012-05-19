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

    for chrom in gdb.get_chromosomes(get_rand=True,
				     get_auto=True,
				     get_sex=True,
				     get_x=True,
				     get_y=True,
				     get_hap=True,
				     get_mito=True):
        print "%s\t%d" % (chrom.name, chrom.length)

main()
