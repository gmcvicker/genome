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
