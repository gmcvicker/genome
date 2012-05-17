
import sys
import numpy as np

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

    sys.stdout.write("%s len:%d n_a:%d n_c:%d n_g:%d n_t:%d n_n:%d "
                     "first_def:%d last_def:%d\n" %
                     (chrom.name, chrom.length, n_a, n_c, n_g, n_t, n_n,
                      node.attrs.first_def_idx + 1,
                      node.attrs.last_def_idx + 1))

    return (n_a, n_c, n_g, n_t)
    
    

def main():
    gdb = genome.db.GenomeDB()

    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <track_name>\n" % sys.argv[0])
        exit(2)

    track_name = sys.argv[1]

    track = gdb.open_track(track_name)

    ttl_a = ttl_c = ttl_t = ttl_g = 0
    
    for chrom in gdb.get_chromosomes():
        (n_a, n_c, n_t, n_g) = get_seq_stats(track, chrom)
        ttl_a += n_a
        ttl_c += n_c
        ttl_g += n_g
        ttl_t += n_t

    sys.stdout.write("total n_a:%d n_c:%d n_g:%d n_t:%d\n" %
                     (ttl_a, ttl_c, ttl_g, ttl_t))
    
    track.close()



if __name__ == "__main__":
    main()
