from collections import Counter
from functools import reduce
from glob_opts import BASES, GAP_PENALTY

import math
import numpy as np


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


def weighted_alphas(raw_counts):
    """Returns the weighted counts for each column. For example, if column
    at position 3 of the alignment only shows some residue (doesn't
    matter which) 150 times, the element of the returned vector at
    index 3 will be 150. The rest of times it shows a gap and is not
    counted.

    """
    return [reduce(lambda res, value: res + value, dic.values(), 0)
            for dic in raw_counts]


def aa_frequencies(counters, nseq):
    """Same structure as raw_counts, but the values are divided by the
    number of sequences.

    """
    return [{residue: count / nseq for residue, count in column.items()}
            for column in counters]


def set_gap_penalties(m, counters):
    """Set the 21st column of the PSSM to the gap penalty that we want for
    that location in the consensus sequence. We set GAP_PENALTY if
    there is little chance of finding a gap and 0 if it is very
    likely.

    """
    for i in range(len(counters)):
        m[20][i] = GAP_PENALTY
        if '-' in counters[i]:
            if counters[i]['-'] > 0.7:
                m[20][i] = - 0.1    # Very little gap penalty, a gap is common

    return m


def pssm_gen(counters, alpha, beta):
    """Returns a PSSM matrix
    """
    swissprotValues = {'A': 8.26, 'R': 5.53, 'N': 4.05, 'D': 5.46, 'C': 1.37,
                       'Q': 3.93, 'E': 6.73, 'G': 7.08, 'H': 2.27, 'I': 5.93,
                       'L': 9.65, 'K': 5.82, 'M': 2.41, 'F': 3.86, 'P': 4.72,
                       'S': 6.60, 'T': 5.35, 'W': 1.09, 'Y': 2.92, 'V': 6.86,
                       'B': 0, 'Z': 0, 'X': 0, '-': 10}

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

    # Get consensus before altering gap values
    consensus = "".join(list(map(lambda x: BASES[np.argmax(x)], np.transpose(m))))
    m = set_gap_penalties(m, counters)

    return {'pssm': m, 'consensus': consensus}
