#!/usr/bin/python

import sys
import tables
import argparse

import genome.db


def parse_args():
    parser = argparse.ArgumentParser(description="write chromosome names "
                                     "and lengths to stdout")

    parser.add_argument("--assembly", default='hg18',
                        help="assembly to print chromosomes for"
                        " (e.g. hg18)")

    parser.add_argument("--rand", action='store_true',
                        help="include random/unknown chromosomes in output",
                        default=False)

    parser.add_argument('--y', action='store_true',
                        help="include chrY in output",
                        default=False)
    
    parser.add_argument('--hap', action='store_true',
                        help="include alternate haplotypes in output",
                        default=False)
        
    parser.add_argument('--mito', action='store_true',
                        help='include mitochondrial chromosomes in output',
                        default=False)

    parser.add_argument('--no_x', action="store_true",
                        default=False,
                        help="do not include chromosome X")
    
    parser.add_argument('--all', action='store_true',
                        help="include all chromosomes in output",
                        default=False)

    parser.add_argument("--ids", action="store_true",
                        default=False,
                        help="print numeric chromosome IDs as well as names")
    
    return parser.parse_args()



def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    
    if args.all:
        chromosomes = gdb.get_all_chromosomes()
    else:
        chromosomes = gdb.get_chromosomes(get_rand=args.rand,
                                          get_y=args.y,
                                          get_x=(not args.no_x),
                                          get_hap=args.hap,
                                          get_mito=args.mito)

    for chrom in chromosomes:
        if args.ids:
            print("%d\t%s\t%d" % (chrom.idnum, chrom.name, chrom.length))
        else:
            print("%s\t%d" % (chrom.name, chrom.length))
        

main()
