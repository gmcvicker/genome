
import sys
import numpy as np

import argparse
import genome.db

import genome.trackstat as trackstat


def set_seq_stats(track, chrom):
    node_name = "/%s" % chrom.name

    if not node_name in track.h5f:
        sys.stderr.write("skipping chromosome %s\n" % chrom.name)
        return

    node = track.h5f.getNode(node_name)

    seq_vals = track.get_nparray(chrom)

    # set counts of each base
    node.attrs.n_a = np.sum(seq_vals == ord('A'))
    node.attrs.n_c = np.sum(seq_vals == ord('C'))
    node.attrs.n_g = np.sum(seq_vals == ord('G'))
    node.attrs.n_t = np.sum(seq_vals == ord('T'))

    # set counts of N and non-N bases
    undef_sites = seq_vals == ord('N')
    node.attrs.n_n = np.sum(undef_sites)
    node.attrs.n_def = seq_vals.size - node.attrs.n_n

    # set index of first and last defined base on chromosome
    w = np.where(seq_vals != ord('N'))[0]

    if w.size == 0:
        raise ValueError("expected at least one defined base on chromosome")

    node.attrs.first_def_idx = w[0]
    node.attrs.last_def_idx = w[-1]

    
    node.flush()
        
    
    
def parse_args():
    parser = argparse.ArgumentParser(description="Computes sequence stats "
                                     "and sets them as attributes on "
                                     "track so that they can be "
                                     "retrieved rapidly (e.g. using "
                                     "get_seq_track_stats.py)")

    parser.add_argument("--assembly", default='hg18',
                        help="assembly to set stats for"
                        " (e.g. hg18)")

    parser.add_argument("--track", default="seq",
                        help="name of sequence track to set stats for")
    

    return parser.parse_args()



def main():    
    args = parse_args()

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    track = gdb.open_track(args.track, "a")

    for chrom in gdb.get_all_chromosomes():
        sys.stderr.write("%s\n" % chrom)
        set_seq_stats(track, chrom)
    
    track.close()



if __name__ == "__main__":
    main()
