
import sys
import os

import numpy as np
import tables
import argparse

import genome.db


UINT16_MAX_VAL = 65535
UINT8_MAX_VAL = 255


def create_carray(track, chrom, dtype):
    if dtype == 'uint8':
        atom = tables.UInt8Atom(dflt=0)
    elif dtype == 'uint16':
        atom = tables.UInt16Atom(dflt=0)
    else:
        raise NotImplementedError("support for dtype %s not "
                                  "yet implemented" % dtype)
    
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray



def combine_tracks(gdb, combined_track, tracks, dtype=None):

    if dtype is None:
        sys.stderr.write("using uint8 datatype by default\n")
        dtype = np.dtype('uint8')

    if dtype == 'uint8':
        max_val = UINT8_MAX_VAL
        larger_dtype = np.uint16
    elif dtype == 'uint16':
        max_val = UINT16_MAX_VAL
        larger_dtype = np.uint32
    else:
        raise NotImplementedError("support for dtype %s not "
                                  "yet implemented" % dtype)
            
    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)
        combined_carray = create_carray(combined_track, chrom, dtype)

        combined = np.zeros(chrom.length, dtype=larger_dtype)

        for track in tracks:
            if track.has_chromosome(chrom):
                combined += track.get_nparray(chrom)

        large_vals = (combined > max_val)

        sys.stderr.write("%d values > max value %d\n" %
                         (np.sum(large_vals), max_val))
        
        combined[large_vals] = max_val
        combined_carray[:] = combined



def create_combined_tracks(combined_track_name, track_names, assembly,
                           dtype=None):
    gdb = genome.db.GenomeDB(assembly=assembly)

    track_list = []
    for track_name in track_names:
        track = gdb.open_track(track_name)
        track_list.append(track)
    
    combined_track = gdb.create_track(combined_track_name)

    combine_tracks(gdb, combined_track, track_list, dtype=dtype)

    for track in track_list:
        track.close()
    
    combined_track.close()




def parse_args():
    parser = argparse.ArgumentParser()
        
    parser.add_argument("--dtype", metavar="", action="store",
                        choices=("uint8", "uint16"), default="uint8",
                        help="datatype of combined track")

    parser.add_argument('--assembly', help="assembly to use", default=None)

    parser.add_argument("combined_track", action="store",
                        help="name of track to store combined counts in")
    
    parser.add_argument("tracks", action="store", nargs="+",
                        help="names of tracks to combine counts from")

    args = parser.parse_args()

    return args



if __name__ == "__main__":
    args = parse_args()
    
    create_combined_tracks(args.combined_track, args.tracks, args.assembly,
                           np.dtype(args.dtype))
        
    
