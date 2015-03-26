import numpy as np
import sys, re, gzip

cimport numpy as np



    
cdef extern from "seq.h":
    ctypedef struct Seq:
        pass
    
    Seq *seq_new()
    long seq_read_fasta_from_file(Seq *seq, const char *filename)
    char *seq_get_seqstr(Seq *seq)
    void seq_free(Seq *seq)

cdef extern from "wig.h":
    float *wig_read_float32(char *filename, long chr_len)
    short *wig_read_int16(char *filename, long chr_len)
    unsigned char *wig_read_uint8(char *filename, long chr_len)

cdef extern from "txtfile.h":
    char *txtfile_read_int8(char *filename, long chr_len,
                            int pos_idx, int val_idx)

    short *txtfile_read_int16(char *filename, long chr_len,
                              int pos_idx, int val_idx)


    
cdef extern from "bedgraph.h":
    short *bedgraph_read_int16(char *filename, long chr_len)
    float *bedgraph_read_float32(char *filename, long chr_len)


cdef extern from "stdlib.h":
    void free(void *ptr)

cdef extern from "string.h":
    void memcpy(void *dest, void *src, size_t n)


cdef extern from "xbf.h":
    ctypedef struct xblist_t:
        pass

    xblist_t *xb_load_mmap(char *filename)
    
    void xb_chrom_vals(xblist_t *xbl, char *chrom_name,
                       unsigned char **fwd_vals,
                       unsigned char **rev_vals, long *chrom_size)

    void xb_free(xblist_t *xbl)







def read_wig_float32(filename, chrom_len):
    cdef np.ndarray[np.float32_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.float32)
    cdef float *vals
    
    vals = wig_read_float32(filename, chrom_len)

    if vals == NULL:
        raise Exception("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(float))
    free(vals)

    return result


def read_wig_int16(filename, chrom_len):
    cdef np.ndarray[np.int16_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.int16)
    cdef short *vals
    
    vals = wig_read_int16(filename, chrom_len)

    if vals == NULL:
        raise Exception("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(short))
    free(vals)

    return result



def read_wig_uint8(filename, chrom_len):
    cdef np.ndarray[np.uint8_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.uint8)
    cdef unsigned char *vals
    
    vals = wig_read_uint8(filename, chrom_len)

    if vals == NULL:
        raise Exception("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(unsigned char))
    free(vals)

    return result




def read_xb(filename, chrom_name, strand):
    cdef np.ndarray[np.uint8_t, ndim=1]result
    cdef unsigned char *fwd_vals
    cdef unsigned char *rev_vals
    cdef long chrom_len
    cdef xblist_t *xbl

    sys.stderr.write("  memory mapping xb file\n")
    xbl = xb_load_mmap(filename)

    sys.stderr.write("  reading chromosome values from xb\n")
    xb_chrom_vals(xbl, chrom_name, &fwd_vals, &rev_vals, &chrom_len)

    result = np.empty(chrom_len, dtype=np.uint8)

    if strand == "forward":
        sys.stderr.write("  copying fwd strand data from xb\n")
        memcpy(result.data, fwd_vals, chrom_len)
    elif strand == "reverse":
        sys.stderr.write("  copying rev strand data from xb\n")
        memcpy(result.data, rev_vals, chrom_len)
    else:
        raise ValueError("unknown strand '%s'" % strand)

    xb_free(xbl)

    return result



def read_bedgraph_int16(filename, chrom_len):
    cdef np.ndarray[np.int16_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.int16)
    cdef short *vals
    
    vals = bedgraph_read_int16(filename, chrom_len)

    if vals == NULL:
        raise IOError("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(short))
    free(vals)

    return result



def read_bedgraph_float32(filename, chrom_len):
    cdef np.ndarray[np.float32_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.float32)
    cdef float *vals
    
    vals = bedgraph_read_float32(filename, chrom_len)

    if vals == NULL:
        raise IOError("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(float))
    free(vals)

    return result



def read_txtfile_int8(filename, chrom_len, pos_idx, val_idx):
    cdef np.ndarray[np.int8_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.int8)
    cdef char *vals

    vals = txtfile_read_int8(filename, chrom_len, pos_idx, val_idx)

    if vals == NULL:
        raise Exception("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(char))
    free(vals)

    return result


def read_txtfile_int16(filename, chrom_len, pos_idx, val_idx):
    cdef np.ndarray[np.int16_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.int16)
    cdef short *vals

    vals = txtfile_read_int16(filename, chrom_len, pos_idx, val_idx)

    if vals == NULL:
        raise Exception("could not read data from file '%s'" % filename)
    
    memcpy(result.data, vals, chrom_len * sizeof(short))
    free(vals)

    return result



def read_fasta(filename, chrom_len):
    cdef np.ndarray[np.uint8_t, ndim=1]result = \
         np.empty(chrom_len, dtype=np.uint8)

    cdef char *c_str
    cdef Seq *seq
    
    seq = seq_new()
    
    sys.stderr.write("reading FASTA from file %s\n" % filename)
    seq_len = seq_read_fasta_from_file(seq, filename)

    if seq_len != chrom_len:
        raise ValueError("expected sequence length to be %d, but "
                         "read %d bp\n" % (chrom_len, seq_len))
    
    c_str = seq_get_seqstr(seq)
    seq_free(seq)

    sys.stderr.write("copying seq memory\n")    
    memcpy(result.data, c_str, chrom_len * sizeof(char))
    free(c_str)

    return result



def read_file(filename, chrom, dtype="float32", format="wiggle",
              pos_idx=None, val_idx=None, strand="forward"):
    """Creates a 1D numpy array of datatype dtype and length
    chrom_len.  Values are read into the array from a file, which
    should be in the specified format (currently 'wiggle', 'bedgraph'
    or 'txtfile').  If the format is 'txtfile' the indices of the
    position and value columns should be provided."""
    
    if format in ("wig", "wiggle"):
        if dtype == "float32":
            return read_wig_float32(filename, chrom.length)
        elif dtype == "int16":
            return read_wig_int16(filename, chrom.length)
        elif dtype == "uint8":
            return read_wig_uint8(filename, chrom.length)
        else:
            raise NotImplementedError("only float32, uint8 and int16 "
                                      "datatypes are currently implemented "
                                      "for wiggle format")
    elif format in ("xb", "xbf"):
        if dtype != "uint8":
            raise NotImplementedError("only uint8 dataype is currently "
                                      "implemented for xb format")
        vals = read_xb(filename, chrom.name, strand)

        if vals.size != chrom.length:
            raise ValueError("length of vector in xb file (%d) "
                             "does not match chromosome length (%d)" %
                             (vals.size, chrom.length))
        return vals
    elif format == "fasta":
        if dtype != "uint8":
            raise NotImplementedError("only uint8 datatype is currently "
                                      "implemented for fasta format")
        
        return read_fasta(filename, chrom.length)
    elif format == "bedgraph":
        if dtype == "int16":
            return read_bedgraph_int16(filename, chrom.length)
        elif dtype == "float32":
            return read_bedgraph_float32(filename, chrom.length)
        else:
            raise NotImplementedError("only int16 and float32 datatypes are "
                                      "currently implemented for bedgraph "
                                      "format")
    elif format == "txtfile":
        if dtype != "int8" and dtype != "int16":
            raise NotImplementedError("only int8 and int16 datatypes are "
                                      "currently implemented for txtfile "
                                      "format")

        if pos_idx is None or pos_idx < 0:
            raise ValueError("pos_idx must be specified in order to read "
                             "txtfile format")
        
        if val_idx is None or val_idx < 0:
            raise ValueError("val_idx must be specified in order to read "
                             "txtfile format")

        if dtype == "int8":
            return read_txtfile_int8(filename, chrom.length, 
                                     pos_idx, val_idx)

        if dtype == "int16":
            return read_txtfile_int16(filename, chrom.length,
                                      pos_idx, val_idx)

        
    else:
        raise NotImplementedError("format '%s' not implemented" % format)
        
