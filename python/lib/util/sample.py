import numpy as np


def uniform_sample(vals, n):
    """Returns an unweighted sample of size n of the provided values"""
    weights = np.ones(len(vals), dtype=np.float32)
    return weighted_sample(vals, weights, n)


def weighted_sample(vals, weights, n):
    """Returns a weighted sample of size n, with replacement, of the
    provided values"""    
    # make cumulative distribution for weight function
    cdf = np.cumsum(weights.astype(np.float32))
    total = cdf[-1]
    cdf = cdf/total
        
    # convert uniform random variates to
    # random variates from weight function, using cumulative distr
    v = np.random.uniform(low=0.0, high=1.0, size=n)
    samp_idx = np.searchsorted(cdf, v)

    # If these values are subsequently used for indexing, it is much
    # faster if the positions are ordered. This may be due to HDF5
    # chunk caching, not sure.
    samp_idx = np.sort(samp_idx)

    # return sampled values
    return vals[samp_idx]
