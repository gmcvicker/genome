#!/usr/bin/python

import sys
import numpy as np
import argparse

import genome.db

import genome.trackstat as trackstat


def get_seq_stats(track, chrom):
    node_name = "/%s" % chrom.name

    if not node_name in track.h5f:
        sys.stderr.write("skipping chromosome %s\n" % chrom.name)
        return

    node = track.h5f.getNode(node_name)

    # set counts of each base
    n_a = node.attrs.n_a
    n_c = node.attrs.n_c
    n_g = node.attrs.n_g
    n_t = node.attrs.n_t

    # set counts of N and non-N bases
    n_n = node.attrs.n_n
    n_def = node.attrs.n_def

    sys.stdout.write("%s len:%d n_a:%d n_c:%d n_g:%d n_t:%d n_n:%d n_def:%d "
                     "first_def:%d last_def:%d\n" %
                     (chrom.name, chrom.length, n_a, n_c, n_g, n_t, n_n,
                      n_def, 
                      node.attrs.first_def_idx + 1,
                      node.attrs.last_def_idx + 1))

    return (n_a, n_c, n_g, n_t, n_def)
    
    

def parse_args():
    parser = argparse.ArgumentParser(description="Write sequence stats "
                                     "to stdout (need to be "
                                     "pre-computed using "
                                     "set_seq_track_stats.py)")

    parser.add_argument("--assembly", default='hg18',
                        help="assembly to get stats for"
                        " (e.g. hg18)")

    parser.add_argument("--track", default="seq",
                        help="name of sequence track to get stats for")
    

    return parser.parse_args()


def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly=args.assembly)

    track = gdb.open_track(args.track)

    ttl_a = ttl_c = ttl_t = ttl_g = ttl_def = 0
    
    for chrom in gdb.get_chromosomes():
        (n_a, n_c, n_t, n_g, n_def) = get_seq_stats(track, chrom)

        ttl_a += n_a
        ttl_c += n_c
        ttl_g += n_g
        ttl_t += n_t
        ttl_def += n_def

    sys.stdout.write("total n_a:%d n_c:%d n_g:%d n_t:%d n_def:%d\n" %
                     (ttl_a, ttl_c, ttl_g, ttl_t, ttl_def))
    
    track.close()



if __name__ == "__main__":
    main()
