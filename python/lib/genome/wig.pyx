import sys

import numpy as np
cimport numpy as np

cdef extern from "wig.h":
     void wig_write_uint8(char *filename, unsigned char *vals, char *chr_name, long chr_len)

ctypedef np.uint8_t DTYPE_t

def write_uint8(filename, np.ndarray[DTYPE_t, ndim=1]vals, chr_name):
     wig_write_uint8(filename, <unsigned char *>vals.data, chr_name, vals.shape[0])

