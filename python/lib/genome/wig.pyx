import sys

import numpy as np
cimport numpy as np

cdef extern from "wig.h":
     void wig_write_uint8(char *filename, unsigned char *vals, char *chr_name, long chr_len)
     void wig_fd_write_uint8(int fd, unsigned char *vals, char *chr_name, long chr_len)
     void wig_write_float32(char *filename, float *vals, char *chr_name, long chr_len)

ctypedef np.uint8_t DTYPE_uint8_t
ctypedef np.float32_t DTYPE_float32_t

def write_fd_uint8(fd, np.ndarray[DTYPE_uint8_t, ndim=1]vals, chr_name):
     wig_fd_write_uint8(fd, <unsigned char *>vals.data, chr_name, vals.shape[0])

def write_uint8(filename, np.ndarray[DTYPE_uint8_t, ndim=1]vals, chr_name):
     wig_write_uint8(filename, <unsigned char *>vals.data, chr_name, vals.shape[0])

def write_float32(filename, np.ndarray[DTYPE_float32_t, ndim=1]vals, chr_name):
    wig_write_float32(filename, <float *>vals.data, chr_name, vals.shape[0])
