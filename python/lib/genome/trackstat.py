
import sys
import numpy as np

import genome.db


class TrackStats(object):
    def __init__(self):
        self.n = 0
        self.n_nan = 0
        self.sum = 0
        self.min = None
        self.max = None


    def mean(self):
        """Calculates mean of sites that are not nan"""
        n = self.n - self.n_nan
        if n == 0:
            return np.inf
        
        return self.sum / float(n)


    def set_from_vals(self, vals):
        self.n = vals.size

        if str(vals.dtype).startswith('float'):
            nan_vals = np.isnan(vals)
            self.n_nan = np.sum(nan_vals)

            if self.n_nan < self.n:
                self.min = np.min(vals[~nan_vals])
                self.max = np.max(vals[~nan_vals])
                self.sum = np.sum(vals[~nan_vals])
        else:
            self.min = np.min(vals)
            self.max = np.max(vals)
            self.sum = np.sum(vals)

        

    def add(self, other):
        self.n += other.n
        self.n_nan += other.n_nan
        self.sum += other.sum
        
        if (self.min is None) or (other.min is not None and 
                                  self.min > other.min):
            self.min = other.min

        if (self.max is None) or (other.max is not None and
                                  self.max < other.max):
            self.max = other.max


    def __str__(self):
        return "n=%d n_nan=%s min=%s max=%s sum=%s" % \
            (self.n, str(self.n_nan), str(self.min), str(self.max), 
             str(self.sum))



def calc_stats(gdb, track):
    """Calculates stats for each chromosome and the entire track,
    but does not store them."""
    
    combined = TrackStats()

    for chrom in gdb.get_chromosomes():
        chrom_stat = TrackStats()        
        vals = track.get_nparray(chrom)
        chrom_stat.set_from_vals(vals)
        sys.stderr.write("%s %s\n" % (str(chrom), str(chrom_stat)))
        combined.add(chrom_stat)

    return combined


def set_stats(gdb, track):
    """Calculates stats for each chromosome and entire track and
    stores them as attributes on the nodes. The provided track must
    be opened in append mode."""
    combined = TrackStats()

    for chrom in gdb.get_all_chromosomes():
        node_name = "/%s" % chrom.name
        if node_name in track.h5f:
            chrom_stat = TrackStats()
            vals = track.get_nparray(chrom)
            chrom_stat.set_from_vals(vals)

            node = track.h5f.getNode("/%s" % chrom.name)
            node.attrs.n = chrom_stat.n
            node.attrs.n_nan = chrom_stat.n_nan
            node.attrs.min = chrom_stat.min
            node.attrs.max = chrom_stat.max
            node.attrs.sum = chrom_stat.sum
            node.flush()

            sys.stderr.write("%s %s\n" % (str(chrom), str(chrom_stat)))
            combined.add(chrom_stat)
        else:
            sys.stderr.write("skipping chromosome %s\n" % chrom)
    
    return combined



def get_stats(gdb, track, chrom=None, verbose=False):
    """Retrieves stats that are stored as attributes. By default
    stats are returned for the whole track, but stats for a
    specific chromosome can also be requested"""
    combined = TrackStats()
    chrom_stat = TrackStats()

    if chrom:
        chrom_list = [chrom]
    else:
        chrom_list = [x for x in gdb.get_chromosomes(get_x=False)]
    
    for chrom in chrom_list:
        node_name = "/%s" % chrom.name
        if node_name in track.h5f:
            node = track.h5f.getNode("/%s" % chrom.name)
            if 'n' not in node.attrs:
                raise ValueError("Stat attributes are not set for track %s"
                                 % track.name)

            chrom_stat.n = node.attrs.n
            chrom_stat.n_nan = node.attrs.n_nan
            chrom_stat.min = node.attrs.min
            chrom_stat.max = node.attrs.max
            chrom_stat.sum = node.attrs.sum

            if verbose:
                sys.stderr.write("%s %s\n" % (str(chrom), str(chrom_stat)))
            combined.add(chrom_stat)

    return combined
        
