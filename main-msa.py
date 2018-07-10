from sequence import *
from pssm import *
from sw import *
from score_matrix import ScoreMatrixSW
from util import msa_parse, seq_parse
from collections import Counter
from functools import reduce

import math
import numpy as np


BASES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

ALIGNMENTS_FILE_NAMES = {"muscle": "data/msaresults-muscle.fasta",
                         "clustal-omega": "data/msaresults-clustalo.fasta",
                         "tcoffee": "data/msaresults-tcoffee.fasta"}

SEQUENCES_FILE_NAME = "data/protein-sequences.fasta"


print("\nMUSCLE PSSM")
alignments = msa_parse(ALIGNMENTS_FILE_NAMES["muscle"])
nseq = len(alignments)
beta = math.sqrt(nseq)
# chain 3 functions like: f(g(h(x))) and obtain the raw_counts
raw_count_per_column = reduce(lambda x, f: f(x), [transpose_aa_chains,
                                                  aa_count,
                                                  remove_gaps], alignments)


alpha = weighted_alphas(raw_count_per_column)
relative_counts = aa_frequencies(raw_count_per_column, nseq)

muscle_pssm_m = pssm(relative_counts, alpha, beta)
# np.set_printoptions(precision=2, threshold=2000) # show full matrix
pssm_in_detail(muscle_pssm_m)


sequences = seq_parse(SEQUENCES_FILE_NAME)
for seq in sequences:
    for pssm in pssms:
        score_matrix = ScoreMatrixSW(seq[1], muscle_pssm_m, -3)
        region = alignSmithWaterman(muscle_pssm_m, score_matrix.s, seq[1], -3)
        print("\nAligning sequence:", seq[0])
        print(region[0])
