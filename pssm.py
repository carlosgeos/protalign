from collections import Counter
from functools import reduce

import math
import numpy as np


BASES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def transpose_aa_chains(alignments):
    """Rows are exchanged for columns. Each element in the list now
    contains all residues that can be found in that column.

    """
    aa_chains = [alignment[1] for alignment in alignments]
    return [''.join(a) for a in zip(*aa_chains)]


def aa_count(columns):
    """The Counter collection automatically counts all the ocurrences of a
    certain element (AA) in 'columns' and stores it in a dict

    """
    return [Counter(column) for column in columns]


def remove_gaps(counters):
    """Ignore '-' in the alignments to perform the frequency calculation

    """
    return [{residue: count for residue, count in column.items() if residue != '-'}
            for column in counters]


def weighted_alphas(raw_counts):
    """Returns the weighted counts for each column. For example, if column
    at position 3 of the alignment only shows some residue (doesn't
    matter which) 150 times, the element of the returned vector at
    index 3 will be 150.

    """
    return [reduce(lambda res, value: res + value, dic.values(), 0)
            for dic in raw_counts]



def aa_frequencies(counters, nseq):
    """Same structure as raw_counts, but the values are divided by the
    number of sequences.

    """
    return [{residue: count / nseq for residue, count in column.items()}
            for column in counters]


def pssm(counters, alpha, beta):
    """Returns a PSSM matrix
    """
    swissprotValues = {'A': 8.26, 'R': 5.53, 'N': 4.05, 'D': 5.46, 'C': 1.37,
                       'Q': 3.93, 'E': 6.73, 'G': 7.08, 'H': 2.27, 'I': 5.93,
                       'L': 9.65, 'K': 5.82, 'M': 2.41, 'F': 3.86, 'P': 4.72,
                       'S': 6.60, 'T': 5.35, 'W': 1.09, 'Y': 2.92, 'V': 6.86,
                       'B': 0, 'Z': 0, 'X': 0}

    nbases = len(BASES)
    ncol = len(counters)
    # q = np.zeros((nBASES, ncol), dtype=float)
    m = np.zeros((nbases, ncol), dtype=float)
    for i in range(ncol):
        for key in counters[i]:
            position = BASES.index(key) # the PSSM has structure
            q = (alpha[i] * counters[i][key] + beta * (swissprotValues[key] / 100)) \
                / (alpha[i] + beta)
            m[position][i] += math.log10(q / (swissprotValues[key] / 100))

    return m


def pssm_in_detail(pssm_m):
    print(pssm_m)
    for i in range(len(pssm_m)):
        print("AminoAcid", BASES[i])
        for j in range(len(pssm_m[0])):
            if pssm_m[i][j] != 0:
                print("Column number", j, end=" : ")
                print(pssm_m[i][j])
