import sys
import tables

import genome.db


def main():
    gdb = genome.db.GenomeDB()

    for chrom in gdb.get_chromosomes():
        print "%s\t%d" % (chrom.name, chrom.length)

main()
