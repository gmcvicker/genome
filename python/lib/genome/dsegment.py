


def find_segments(scores, thresh, max_dropoff):
    """Returns a set of high-scoring segments for the provided sequence
    of scores. The segments returned are half-open (start, end) tuples,
    corresponding to the indices of the provided scores. The returned
    segments have scores of at least thresh with no drop greater than
    max_dropoff within the segment."""
    start = end = max_score = cur_score = i = 0
    segments = []

    for score in scores:
        cur_score += score
        if cur_score >= max_score:
            max_score = cur_score
            end = i+1

        if (cur_score < 0) or (cur_score < max_score - max_dropoff):
            if max_score >= thresh:
                # create a new segment
                seg = (start,end)
                segments.append(seg)

            max_score = cur_score = 0
            start = end = i+1
        
        i += 1

    # handle final segment
    if max_score >= thresh:
        seg = (start, end)
        segments.append(seg)

    return segments

