from sequence import Sequence
from pssm import transpose_aa_chains, aa_count, \
    weighted_alphas, aa_frequencies, pssm_gen
from sw_msa import align_smith_waterman
from backtrack_matrix_msa import BacktrackMatrixSWMSA
from util import msa_parse, seq_parse
from functools import reduce
from glob_opts import ALIGNMENTS_FILE_NAMES, FULL_SEQUENCES_FILE, BASES, GAP_PENALTY

import math
import numpy as np



pssm_collection = {}
sequences = [Sequence(seq) for seq in seq_parse("data/" + FULL_SEQUENCES_FILE + ".fasta")]

for name, file_path in ALIGNMENTS_FILE_NAMES.items():
    alignments = msa_parse(ALIGNMENTS_FILE_NAMES[name])
    nseq = len(alignments)
    beta = math.sqrt(nseq)
    # chain 2 functions like: f(g(x)) and obtain the raw_counts
    raw_count_per_column = reduce(lambda x, f: f(x), [transpose_aa_chains,
                                                      aa_count], alignments)

    alpha = weighted_alphas(raw_count_per_column)
    relative_counts = aa_frequencies(raw_count_per_column, nseq)

    pssm_info = pssm_gen(relative_counts, alpha, beta)
    pssm = np.transpose(pssm_info['pssm'])
    consensus = pssm_info['consensus']
    pssm_collection[name] = pssm

    print(name + " PSSM\n")
    print(pssm)

    print("Consensus of PSSM", name, ":")
    print(consensus, end="\n\n\n")

    for seq in sequences:
        backtrack_matrix = BacktrackMatrixSWMSA(seq, pssm)
        print("\n\n\n\n\n")
        # print(backtrack_matrix.s[267][34])
        # for row in backtrack_matrix.s:
        #     for elem in row:

        align_smith_waterman(seq, backtrack_matrix.s, consensus)
