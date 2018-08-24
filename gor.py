from glob_opts import *
from math import log

def fill_counters(sequences, counters):
    """Fill the three-dimensional array with observations of aminoacids
    and conformations.

    It contains positional information according to the indices in
    CONFORMATIONS and BASES

    """
    f_sr_i = len(BASES)         # index to store f_sr
    for seq in sequences:
        for i in range(len(seq)):
            aa_index = BASES.index(seq["aa", i])
            conf_index = CONFORMATIONS.index(seq["conf", i])
            counters[conf_index][aa_index][f_sr_i] += 1
            for m in seq.neighbours(i):
                counters[conf_index][aa_index][m] += 1 # f_srmr
    return counters


def self_information(s_index, r_index, counters):
    """Gives the self information between S and R

    """

    ns = [x for x in CONFORMATIONS_INDICES if x != s_index] # not s !

    f_sr = counters[s_index][r_index][f_sr_i]
    f_nsr = counters[ns[0]][r_index][f_sr_i] + \
        counters[ns[1]][r_index][f_sr_i]

    f_ns = 0
    for r in counters[ns[0]]:
        f_ns += r[f_sr_i]
    for r in counters[ns[1]]:
        f_ns += r[f_sr_i]
    f_s = 0
    for r in counters[s_index]:
        f_s += r[f_sr_i]

    return log(f_sr / f_nsr) + log(f_ns / f_s)


def cond_information(s_index, r_index, rm_index, counters):
    """Calculates the conditional information for a given conformation and
    residue, taking into account its neighbourhood

    """
    ns = [x for x in CONFORMATIONS_INDICES if x != s_index] # not s !

    f_srmr = counters[s_index][r_index][rm_index]
    f_nsrmr = counters[ns[0]][r_index][rm_index] + \
        counters[ns[1]][r_index][rm_index]
    f_nsr = counters[ns[0]][r_index][f_sr_i] + \
        counters[ns[1]][r_index][f_sr_i]
    f_sr = counters[s_index][r_index][f_sr_i]

    first_operand = float("-inf")
    if f_srmr != 0 and f_nsrmr != 0:
        first_operand = log(f_srmr / f_nsrmr)

    second_operand = float("-inf")
    if f_nsr != 0 and f_sr != 0:
        second_operand = log(f_nsr / f_sr)

    return first_operand + second_operand


def total_information(s_index, r_index, neighbours, counters):
    """Computes the total amount of information for conformation S and
    residue R as described in the GOR III algorithm

    """
    self_info = self_information(s_index, r_index, counters)
    total_pair_information = 0
    for n in neighbours:
        total_pair_information += cond_information(s_index, r_index, n, counters)

    return self_info + total_pair_information


def gor(seq, counters):
    """Implements GOR and returns a list of predictions given a list of
    sequences and pre-calculated counters (training)

    """
    pred = ""
    for i in range(len(seq)):
        conf_scores = [0 for conf in CONFORMATIONS]
        neighbours = seq.neighbours(i)
        r_index = BASES.index(seq["aa", i])
        for s_index in range(len(CONFORMATIONS)):
            conf_scores[s_index] = total_information(s_index, r_index,
                                                         neighbours, counters)
        the_best = 0
        for i in range(1, len(conf_scores)):
            if conf_scores[i] > conf_scores[the_best]:
                the_best = i
        pred += CONFORMATIONS[the_best]
    return pred
