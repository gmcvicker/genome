import sys

import numpy as np
cimport numpy as np



def find_segments(np.ndarray[np.int16_t, ndim=1] scores,
                  int thresh, int max_dropoff):
    """Returns a set of high-scoring segments for the provided sequence
    of scores. The segments returned are half-open (start, end) tuples,
    corresponding to the indices of the provided scores. The returned
    segments have scores of at least thresh with no drop greater than
    max_dropoff within the segment."""
    cdef int start, end, max_score, cur_score, i, chr_len

    i = 0
    cur_score = 0
    max_score = 0
    end = 0
    start = 0

    chr_len = len(scores)
    
    segments = []

    for i in xrange(chr_len):
        if (i % 1000000) == 0:
            sys.stderr.write("." )

        cur_score += scores[i]
        if cur_score >= max_score:
            max_score = cur_score
            end = i+1

        if (cur_score < 0) or (cur_score < max_score - max_dropoff):
            if max_score >= thresh:
                # create a new segment
                seg = (start,end)
                segments.append(seg)
                #sys.stderr.write("segment: %d %d %d %d\n" %
                #                 (start+1, end, max_score, end-start))

            max_score = cur_score = 0
            start = end = i+1
        
        i += 1

    # handle final segment
    if max_score >= thresh:
        seg = (start, end)
        segments.append(seg)

    sys.stderr.write("\n")
        
    return segments

