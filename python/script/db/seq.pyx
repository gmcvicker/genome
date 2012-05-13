import sys

cdef extern from "zlib.h":
    ctypedef void * gzFile
    gzFile gzopen(char *filename, char *mode)
    gzclose(gzFile f)


cdef extern from "seq.h":
    struct Seq:
        pass

    Seq *seq_new()
    int seq_read_fasta_record(Seq *seq, gzFile f)
    char *seq_get_seqstr(Seq *seq)
    void seq_free(Seq *seq)


def read_fasta(filename):
    cdef Seq *seq
    cdef gzFile gzf
    
    gzf = gzopen(filename, "rb")
    seq = seq_new()
    seq_read_fasta_record(seq, gzf)
    seq_str = seq_get_seqstr(seq)
    seq_free(seq)

    # gzclose(gzf)

    return seq_str
