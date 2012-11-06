
import sys
import os
import numpy as np

import genome.db
import genome.seq
cimport numpy as np



def get_snp_kmers(np.ndarray[np.uint8_t, ndim=1] kmer_array,  
                  np.ndarray[np.int64_t, ndim=1] snp_idx, 
                  np.ndarray[np.uint8_t, ndim=1] ref_base, 
                  np.ndarray[np.uint8_t, ndim=1] alt_base, 
                  int i=0, kmer_list=[]):

    idx = snp_idx[i]

    # append kmer with reference snp at this position
    kmer_array[idx] = chr(ref_base[idx])
    kmer_list.append("".join(kmer_array))
    
    if i < snp_idx.size-1:
        # recursively add other SNPs
        get_snp_kmers(kmer_array, snp_idx, ref_base, alt_base, i=i+1,
                      kmer_list=kmer_list)

    
    # append kmer with alt snp at this position
    kmer_array[idx] = chr(alt_base[idx])
    kmer_list.append("".join(kmer_array))

    if i < snp_idx.size-1:
        # recursively add other SNPs
        get_snp_kmers(kmer_array, snp_idx, ref_base, alt_base, i=i+1,
                      kmer_list=kmer_list)

    return kmer_list

    


def count_kmers_with_snps(dna, int kmer_size, kmer_counts, 
                          np.ndarray[np.uint8_t, ndim=1] snp_ref,
                          np.ndarray[np.uint8_t, ndim=1] snp_alt):
    cdef int i, chr_len, n_obs
    chr_len = len(dna)

    for i in xrange(chr_len-kmer_size+1):
        if (i % 1000000) == 0:
            sys.stderr.write(".")

        kmer = dna[i:(i + kmer_size)]

        # account for SNPs here!
        ref = snp_ref[1:(i + kmer_size)]
        alt = snp_ref[1:(i + kmer_size)]

        snp_idx = np.where(ref != 0)[0]
        

        kmer_list = []
        if snp_idx.size > 0:
            # there are SNPs that overlap this kmer
            kmer_array = np.array([x for x in kmer])
            
            kmer_list = get_snp_kmers(kmer_array, snp_idx, ref, alt)

            sys.stderr.write("kmer overlaps SNP!:\n  " +
                             "\n  ".join(kmer_list))
            
        else:
            kmer_list = [kmer]

        # update the count of kmers from all sites

        for k in kmer_list:
            if kmer in kmer_counts:
                kmer_counts[k] += 1
            else:
                kmer_counts[k] = 1
            
    return





def count_kmers(dna, int kmer_size, np.ndarray[np.int32_t, ndim=1] obs_counts,
                all_kmers, obs_kmers):
    """Simultneously Counts kmers at all sites and at a subset of
    flagged sites with counts"""
    cdef int i, chr_len, n_obs
    chr_len = obs_counts.shape[0]

    for i in xrange(chr_len-kmer_size+1):
        if (i % 1000000) == 0:
            sys.stderr.write(".")

        # update the count of kmers from all sites
        kmer = dna[i:(i + kmer_size)]
        if kmer in all_kmers:
            all_kmers[kmer] += 1
        else:
            all_kmers[kmer] = 1

        # is this site flagged with a count? if so count it as 'observed'
        n_obs = obs_counts[i]
        if n_obs > 0:
            # update the count of kmers from observed cuts
            if kmer in obs_kmers:
                obs_kmers[kmer] += n_obs
            else:
                obs_kmers[kmer] = n_obs
    return



def count_observed_kmers(dna, int kmer_size, int chr_len,
                         np.ndarray[np.int32_t, ndim=1] kmer_starts,
                         kmer_counts):

    cdef int start
    
    for start in kmer_starts:
        end = start + kmer_size - 1
        
        if end > chr_len:
            # off end of chromosome
            continue
            
        kmer = dna[(start-1):end]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1

    return
        


def count_observed_joint_kmers(dna, int kmer_size, int kmer_offset,
                               np.ndarray[np.int32_t, ndim=1] frag_start,
                               np.ndarray[np.int32_t, ndim=1] frag_size,
                               obs_kmers):
    cdef int i, n_frag, chr_len, left_start, left_end, right_start, right_end

    n_frag = frag_start.shape[0]
    chr_len = len(dna)

    for i in xrange(n_frag):
        # get sequence of kmer at left cutsite
        left_start = frag_start[i] + kmer_offset - 1
        left_end = left_start + kmer_size
        if left_start >= 0 and left_end <= chr_len:
            left_kmer = dna[left_start:left_end]
        else:
            continue

        # get sequence of kmer at right cutsite
        right_end = frag_start[i] + frag_size[i] - kmer_offset - 1
        right_start = right_end - kmer_size
        if right_start >= 0 and right_end <= chr_len:
            right_kmer = dna[right_start:right_end]
        else:
            continue

        key = left_kmer + ":" + right_kmer

        if key in obs_kmers:
            obs_kmers[key] += 1
        else:
            obs_kmers[key] = 1

        if i % 100000 == 0:
            sys.stderr.write("  sample: %d, frag_start: %d, frag_size: %d\n" %
                             (i, frag_start[i], frag_size[i]))


    return





    
def count_joint_kmers(dna, int kmer_size, int kmer_offset, int kmer_sep,
                      np.ndarray[np.uint8_t, ndim=1] obs_fwd_counts,
                      all_kmers, obs_kmers):
    
    """Counts joint occurances of kmers that are separated by a fixed
    difference.  kmer_offset - offset of the start of the kmer from
    the 'observed' positions.  kmer_sep - separation between the
    forward observed and reverse positions.  obs_fwd_counts - array
    (length of chromosome), with counts of forward observations
    """    
    cdef int i, chr_len, n_obs, start, end, min_idx, max_idx, count
    
    chr_len = obs_fwd_counts.shape[0]

    min_idx = min(0, -kmer_offset)
    max_idx = chr_len + kmer_offset - kmer_sep

    sys.stderr.write("scanning between %d - %d\n" % (min_idx+1, max_idx+1))

    count = 0
    for i in xrange(min_idx, max_idx):
        start =  i + kmer_offset
        end = start + kmer_size
        left_kmer = dna[start:end]

        count += 1
        if count > 1000000:
            sys.stderr.write(".")
            count = 0

        end   = i + kmer_sep - kmer_offset + 1
        start = end - kmer_size
        right_kmer = dna[start:end]

        key = left_kmer + ":" + right_kmer

        if key in all_kmers:
            all_kmers[key] += 1
        else:
            all_kmers[key] = 1

        if obs_fwd_counts[i]:
            if key in obs_kmers:
                obs_kmers[key] += obs_fwd_counts[i]
            else:
                obs_kmers[key] = obs_fwd_counts[i]

    return
        
    



def list_kmers(kmer_size, prefix_list=[""]):
    """Returns a list of all possible DNA kmers of a given size"""
    kmer_list = []

    if kmer_size < 1:
        return prefix_list
    
    new_pfx_list = []
    for pfx in prefix_list:
        for nuc in ('A', 'C', 'G', 'T'):
            new_pfx_list.append(nuc + pfx)

    return list_kmers(kmer_size-1, prefix_list=new_pfx_list)



def assign_rates(dna, int kmer_size, int kmer_offset, kmer_cut_rates,
                 np.ndarray[np.float32_t, ndim=1] rate):
    """Assigns expected cutting rates at each sequence position, based on
    a dictionary of rate kmers"""
    cdef int i, chr_len, n_obs
    chr_len = rate.shape[0]

    for i in xrange(chr_len-kmer_size+1):
        if (i % 1000000) == 0:
            sys.stderr.write(".")

        kmer = dna[i:(i + kmer_size)]

        if kmer in kmer_cut_rates:
            rate[i + kmer_offset] = kmer_cut_rates[kmer]
        else:
            rate[i + kmer_offset] = np.nan

    return



def set_joint_mnase_weights(dna, kmer_rate_dict, kmer_size, kmer_offset,
                            np.ndarray[np.float32_t, ndim=1] size_distr,
                            np.ndarray[np.float32_t, ndim=1] mnase_expect,
                            int min_frag_size, int max_frag_size):
    cdef int chrom_len, i, j, n_size, is_odd
    cdef int left_cut_idx, right_cut_idx, frag_size, offset
    cdef float size_prob, rate

    # range of MNase fragment sizes we are considering
    n_size = max_frag_size - min_frag_size + 1

    chrom_len = len(dna)

    # loop over every base on chromosome
    for i in xrange(chrom_len):
        # loop over every possible fragment size centered at this base

        if (i % 1000000) == 0:
            sys.stderr.write(".")

        # set flag to true if first fragment size is odd
        is_odd = (min_frag_size % 2) > 0
        ttl_rate = 0.0

        for j in xrange(n_size):
            frag_size = min_frag_size + j

            # weight each fragment size using size distribution
            size_prob = size_distr[j]

            if is_odd:
                offset = (frag_size - 1) / 2
                left_start = i - offset + kmer_offset
                left_end = left_start + kmer_size

                right_end = i + offset - kmer_offset + 1
                right_start = right_end - kmer_size

                if (left_start >= 0) and (right_end <= chrom_len):
                    left_kmer = dna[left_start:left_end]
                    right_kmer = dna[right_start:right_end]

                    kmer_key = left_kmer + ":" + right_kmer

                    if kmer_key in kmer_rate_dict:
                        ttl_rate += kmer_rate_dict[kmer_key] * size_prob
            else:
                offset = frag_size / 2

                # For even fragment sizes, there are two "half" midpoints instead
                # of one. Consider both of them, but weight by 0.5.
                left_start = i - offset + kmer_offset
                left_end   = left_start + kmer_size
                right_end = i + offset - kmer_offset + 1
                right_start = right_end - kmer_size
                
                if (left_start >= 0) and (right_end <= chrom_len):
                    left_kmer = dna[left_start:left_end]
                    right_kmer = dna[right_start:right_end]
                    kmer_key = left_kmer + ":" + right_kmer
                    if kmer_key in kmer_rate_dict:
                        # add half the rate
                        ttl_rate += 0.5 * kmer_rate_dict[kmer_key] * size_prob

                # shift one over for the other fragment we are considering
                left_start += 1
                left_end   += 1
                right_end  += 1
                right_start += 1
                
                if right_end <= chrom_len:
                    left_kmer = dna[left_start:left_end]
                    right_kmer = dna[right_start:right_end]
                    kmer_key = left_kmer + ":" + right_kmer
                    if kmer_key in kmer_rate_dict:
                        # add another half rate for the other possible fragment
                        ttl_rate += 0.5 * kmer_rate_dict[kmer_key] * size_prob
                
            # alternate between odd and even fragment sizes
            is_odd = 0 if is_odd else 1

        mnase_expect[i] = ttl_rate




def set_mnase_weights(np.ndarray[np.float32_t, ndim=1] fwd_expect_cuts,
                      np.ndarray[np.float32_t, ndim=1] rev_expect_cuts,
                      np.ndarray[np.float32_t, ndim=1] size_distr,
                      np.ndarray[np.float32_t, ndim=1] mnase_expect,
                      int min_size, int max_size):
    cdef int chrom_len, i, j, n_size, is_odd
    cdef int left_cut_idx, right_cut_idx, size, offset
    cdef float size_prob, left_cut_rate, right_cut_rate, ttl_expect_rate

    # range of MNase fragment sizes we are considering
    n_size = max_size - min_size + 1

    chrom_len = fwd_expect_cuts.shape[0]

    # loop over every base on chromosome
    for i in xrange(chrom_len):
        # loop over every possible fragment size centered at this base
        ttl_expect_rate = 0.0

        if (i % 1000000) == 0:
            sys.stderr.write(".")

        # set flag to true if first fragment size is odd
        is_odd = (min_size % 2) > 0
        for j in xrange(n_size):
            size = min_size + j

            # weight each fragment size using size distribution
            size_prob = size_distr[j]

            if is_odd:
                offset = (size - 1) / 2
                # expected cut rates are at base immediately 5' of cutsite
                left_cut_idx = i - offset - 1
                right_cut_idx = i + offset + 1

                if (left_cut_idx >= 0) and (right_cut_idx < chrom_len):
                    left_cut_rate = fwd_expect_cuts[left_cut_idx]
                    right_cut_rate = rev_expect_cuts[right_cut_idx]

                    # How to get expected rate? Several possible options, but not
                    # sure what is correct. For now assume that we can just use mean
                    # of rates from each end of fragment
                    ttl_expect_rate += size_prob * 0.5 * (left_cut_rate + right_cut_rate)
            else:
                offset = size / 2

                # For even fragment sizes, there are two "half" midpoints instead
                # of one. Consider both of them, but weight by 0.5.
                left_cut_idx = i - offset - 1
                right_cut_idx = i + offset
                if (left_cut_idx >= 0) and (right_cut_idx < chrom_len):
                    left_cut_rate = fwd_expect_cuts[left_cut_idx]
                    right_cut_rate = rev_expect_cuts[right_cut_idx]
                    ttl_expect_rate += size_prob * 0.25 * (left_cut_rate + right_cut_rate)

                left_cut_idx = i - offset
                right_cut_idx = i + offset + 1
                if (left_cut_idx >= 0) and (right_cut_idx < chrom_len):
                    left_cut_rate = fwd_expect_cuts[left_cut_idx]
                    right_cut_rate = rev_expect_cuts[right_cut_idx]
                    ttl_expect_rate += size_prob * 0.25 * (left_cut_rate + right_cut_rate)

            # alternate between odd and even fragment sizes
            is_odd = 0 if is_odd else 1

        mnase_expect[i] = ttl_expect_rate
