import sys
import tables

import genome.db


def main():
    gdb = genome.db.GenomeDB()

    for chrom in gdb.get_chromosomes(get_rand=True,
				     get_auto=True,
				     get_sex=True,
				     get_x=True,
				     get_y=True,
				     get_hap=True,
				     get_mito=True):
        print "%s\t%d" % (chrom.name, chrom.length)

main()
