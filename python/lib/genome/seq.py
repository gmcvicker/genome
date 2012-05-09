
import string
import numpy as np

dna_comp = None

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global dna_comp

    if dna_comp is None:
        dna_comp = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(dna_comp)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]


def revcomp_nparray(vals):
    seqstr = from_nparray(vals)
    seqstr = revcomp(seqstr)
    return np.array([ord(x) for x in seqstr], dtype=np.uint8)
    

def from_nparray(vals):
    """converts a numpy array into a sequence string"""
    return "".join(chr(x) for x in vals)
