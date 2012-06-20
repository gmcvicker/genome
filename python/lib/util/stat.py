import numpy as np


def assign_quantiles(val_list, n_quantile=5):
    vals = np.array(val_list, dtype=np.float32)
    vals.sort()

    q = float(vals.size) / float(n_quantile)

    quantiles = []

    # compute value of each quantile
    for i in range(n_quantile-1):
        idx = int(round(q * float(i+1)))

        if idx >= vals.size:
            idx = vals.size - 1

        quantiles.append(vals[idx])

    sys.stderr.write("quantiles:\n"
                     "  %s\n" % (", ".join([str(x) for x in quantiles])))

    # assign values to each quantile
    q_vals = np.zeros(vals.size, dtype=np.int8)
    for i in range(len(val_list)):
        for j in range(len(quantiles)):
            if val_list[i] < quantiles[j]:
                break

        if val_list[i] > quantiles[j]:
            j += 1
        q_vals[i] = j
    
    return q_vals



def split_sample(counts, frac=0.5):
    """Returns a subsample of the counts at each position in the provided
    array of counts"""
    n = np.sum(counts)
    index = np.zeros(n, dtype=np.uint16)

    # put indexes into array n times, where n is count at index
    j = 0
    for i in np.where(counts)[0]:
        count = counts[i]
        index[j:j+count] = i
        j += count
    
    # shuffle the indexes
    np.random.shuffle(index)
    half = int(n * frac)


    if half < 1:
        # there are very few datapoints in this region
        return np.zeros(len(counts), dtype=counts.dtype)

    # count the number of times each index occurs in the first
    # frac of the shuffled set. This is our subsample.
    return np.bincount(index[0:half], minlength=len(counts))


        
        
def match_samples(samp1, samp2):
    """Downsamples the larger of two samples so that they match"""
    n1 = np.sum(samp1)
    n2 = np.sum(samp2)

    if n1 > n2:
        frac = float(n2) / float(n1)
        samp1 = split_sample(samp1, frac)
    elif n2 > n1:
        frac = float(n1) / float(n2)
        samp2 = split_sample(samp2, frac)

    return (samp1, samp2)
