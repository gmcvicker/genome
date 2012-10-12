import numpy as np

NUCS = ('A', 'C', 'G', 'T')

NUC_ID_UNDEF = -1
NUC_ID_A = 0
NUC_ID_C = 1
NUC_ID_G = 2
NUC_ID_T = 3
NUC_ID_N = -1
N_NUC = 4

DINUC_ID_UNDEF = -1
N_DINUC = 16

AA_ID = 0
TT_ID = 1
NON_AA_TT_ID = 2
N_AA_TT = 3


def build_dinuc_vector(dna_seq):
    return build_dinuc_matrix(dna_seq).flatten()


def build_dinuc_matrix(dna_seq):
    matrix = np.zeros((len(dna_seq)-1, N_DINUC), dtype=np.uint8)

    for i in range(len(dna_seq)-1):
        dinuc_id = dinuc2id(dna_seq[i], dna_seq[i+1])

        if dinuc_id == NUC_ID_UNDEF:
            return None

        matrix[i, dinuc_id] = 1
    
    return matrix


def build_aa_tt_matrix(dna_seq):
    matrix = np.zeros((len(dna_seq)-1, N_AA_TT), dtype=np.uint8)

    for i in range(len(dna_seq)-1):
        if dna_seq[i] == "A" and dna_seq[i+1] == "A":
            # this is an AA dinucleotide
            matrix[i, AA_ID] = 1
        elif dna_seq[i] == "T" and dna_seq[i+1] == "T":
            # this is a TT dinucleotide
            matrix[i, TT_ID] = 1
        else:
            # this is neither an AA nor a TT dinucleotide
            matrix[i, NON_AA_TT_ID] = 1

    return matrix


def correct_aa_tt_matrix(aa_tt_matrix):
    mean_vals = np.zeros(N_NUC, dtype=np.float32)

    for i in range(N_AA_TT):
        mean_vals = np.mean(aa_tt_matrix[:,i])

    return aa_tt_matrix - mean_vals


def correct_nuc_matrix(nuc_matrix):
    """Corrects a nucleotide matrix, by subtracting a vector of
    mean nucleotide proportions from each row"""

    nuc_mean = np.zeros(N_NUC, dtype=np.float32)
    for i in range(N_NUC):
        nuc_mean[i] = np.mean(nuc_matrix[:,i])

    return nuc_matrix - nuc_mean
    

def correct_dinuc_matrix(nuc_matrix, dinuc_matrix):
    # calculate mono-nucleotide frequencies for this sequence
    nuc_mean = np.zeros(N_NUC, dtype=np.float32)
    for i in range(N_NUC):
        nuc_mean[i] = np.mean(nuc_matrix[:,i])

    # now calculate expected dinucleotide frequencies, assuming independence of
    # mono-nucleotides (this indepdendence assumption is wrong, but hopefully
    # good enough for this correction)
    dinuc_expect = np.zeros(N_DINUC, dtype=np.float32)
    for i in range(N_NUC):
        for j in range(N_NUC):
            dinuc_id = (i * N_NUC) + j
            dinuc_expect[dinuc_id] = nuc_mean[i] * nuc_mean[j]

    return dinuc_matrix - dinuc_expect

    

def build_nuc_vector(dna_seq):
    return build_nuc_matrix(dna_seq).flatten()



def build_nuc_matrix(dna_seq):
    matrix = np.zeros((len(dna_seq), N_NUC), dtype=np.uint8)

    for i in range(len(dna_seq)):
        nuc_id = nuc2id(dna_seq[i])

        if nuc_id == NUC_ID_UNDEF:
            return None

        matrix[i, nuc_id] = 1

    return matrix



def nuc2id(nuc):
    if nuc == "A" or nuc == "a":
        return NUC_ID_A
    if nuc == "C" or nuc == "c":
        return NUC_ID_C
    if nuc == "G" or nuc == "g":
        return NUC_ID_G
    if nuc == "T" or nuc == "t":
        return NUC_ID_T

    return NUC_ID_UNDEF


def id2nuc(nuc_id):
    if nuc_id == NUC_ID_A:
        return "A"
    if nuc_id == NUC_ID_C:
        return "C"
    if nuc_id == NUC_ID_G:
        return "G"
    if nuc_id == NUC_ID_T:
        return "T"

    return "N"



def dinuc2id(nuc1, nuc2):
    id1 = nuc2id(nuc1)
    id2 = nuc2id(nuc2)

    if id1 == NUC_ID_UNDEF:
        return NUC_ID_UNDEF
    if id2 == NUC_ID_UNDEF:
        return NUC_ID_UNDEF

    return (id1 * N_NUC) + id2


def id2dinuc(dinuc_id):
    nuc_id2 = dinuc_id % N_NUC
    nuc_id1 = (dinuc_id - nuc_id2) / N_NUC

    return (id2nuc(nuc_id1), id2nuc(nuc_id2))



def get_all_nuc():
    return ['A', 'C', 'T', 'G']

def get_all_dinuc():
    dinuc = []
    for i in range(N_DINUC):
        nuc1, nuc2 = id2dinuc(i)
        dinuc.append(nuc1 + nuc2)
    return dinuc
