import sys

import numpy as np
cimport numpy as np

# use minimum int32 value to indicate an undefined distance
cdef int UNDEF_DIST = -2147483648


def calc_dists_tuples_with_strand(start_end_strand_list, chr_len):
    """Returns a 1D numpy array of int32s that is the length of the
    chromosome. The array contains distances to the closest feature for every
    position on the chromosome (element 0 of the array is chromosome
    position 1). Features must be provided as list of (start,end,strand)
    tuples.  Distances are signed by whether the nearest feature is
    'upstream' or 'downstream' to that position. Positions where a
    downstream feature is closest are given negative distances; positions
    where an upstream feature is closest are given positive distances.
    A feature is considered downstream of a position if it is to the right
    of the position and on the forward strand OR if it is to the left of
    the position and on the reverse strand."""
    fwd_dists = np.empty(chr_len, dtype=np.int32)
    rev_dists = np.empty(chr_len, dtype=np.int32)

    fwd_dists[:] = UNDEF_DIST
    rev_dists[:] = UNDEF_DIST

    for x in start_end_strand_list:
        if x[2] == 1:
            fwd_dists[(x[0]-1):x[1]] = 0
        elif x[2] == -1:
            rev_dists[(x[0]-1):x[1]] = 0
        else:
            raise ValueError("unknown strand '%d'" % x[2])

    calc_dists_array(fwd_dists, flag_direction=True)
    calc_dists_array(rev_dists, flag_direction=True)

    # flip sign of distances for reverse strand elements
    # because up/downstream is flipped
    rev_dists = -rev_dists

    # take whichever absolute distance is lower from fwd/rev elements
    rev_dist_lower = np.abs(rev_dists) < np.abs(fwd_dists)
    fwd_dists[rev_dist_lower] = rev_dists[rev_dist_lower]

    return fwd_dists



def calc_dists_tuples(start_end_list, chr_len,
                      flag_direction=False):
    """Given a set of features represented by a list of (start,end)
    tuples, returns the physical distances along the chrosomes to
    each of the features. The distances are returned as 1D numpy array
    of int32s that is the length of the chromosome."""
    dists = np.empty(chr_len, dtype=np.int32)

    # set all distances to UNDEFINED
    dists[:] = UNDEF_DIST
    for start_end in start_end_list:
        # flag regions defined by start_end tuples as distance 0
        dists[(start_end[0]-1):start_end[1]] = 0

    calc_dists_array(dists, flag_direction=flag_direction)
    
    return dists




def calc_dists_array(np.ndarray[np.int32_t, ndim=1] dists,
                     flag_direction=False):
    """Calculates physical distances to features along a chromosome.
    The distances are set in the provided numpy array, which must be
    initialized so that the locations of the features are flagged with
    0 and all other elements have value UNDEF_DIST. If flag_direction is True
    then distances to the left of features are given negative values, and
    distances to the right of features are given positive values. The
    shortest distance is still considered to be the smallest absolute
    value."""
    cdef int i, rev_i, chr_len, cur_dist

    cur_dist = UNDEF_DIST
    chr_len = dists.shape[0]
    
    for i in xrange(chr_len):
        if dists[i] == 0:
            cur_dist = 0
        elif cur_dist >= 0:
            cur_dist += 1
        dists[i] = cur_dist

    # move backwards on chrom, calcing distance from right coords
    cur_dist = 1
    for i in xrange(chr_len):
        rev_i = chr_len - i - 1

        if dists[rev_i] == 0:
            cur_dist = 0
        elif cur_dist <= 0:
            cur_dist -= 1

            if (dists[rev_i] == UNDEF_DIST) or (abs(cur_dist) < dists[rev_i]):
                # update, distance from right element is shorter
                dists[rev_i] = cur_dist

    if not flag_direction:
        dists[:] = np.abs(dists)

    return dists
