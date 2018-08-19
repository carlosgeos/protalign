from sequence import Sequence
from pssm import transpose_aa_chains, aa_count, \
    weighted_alphas, aa_frequencies, pssm_gen
from sw_msa import align_smith_waterman
from backtrack_matrix_msa import BacktrackMatrixSWMSA
from util import msa_parse, seq_parse
from functools import reduce
from glob_opts import ALIGNMENTS_FILE_NAMES, FULL_SEQUENCES_FILE, BASES

import math
import numpy as np



pssm_collection = {}
for name, file_path in ALIGNMENTS_FILE_NAMES.items():
    alignments = msa_parse(ALIGNMENTS_FILE_NAMES[name])
    nseq = len(alignments)
    beta = math.sqrt(nseq)
    # chain 3 functions like: f(g(h(x))) and obtain the raw_counts
    raw_count_per_column = reduce(lambda x, f: f(x), [transpose_aa_chains,
                                                      aa_count], alignments)

    alpha = weighted_alphas(raw_count_per_column)
    relative_counts = aa_frequencies(raw_count_per_column, nseq)
    print(relative_counts[40])

    pssm = np.transpose(pssm_gen(relative_counts, alpha, beta))
    pssm_collection[name] = pssm

    print(name + " PSSM\n")
    print(pssm)
    print(len(pssm))

    consensus = "".join(list(map(lambda x: BASES[np.argmax(x)], pssm)))
    print("Consensus of PSSM", name, ":")
    print(consensus, end="\n\n\n")

test = pssm_collection["muscle"]

print(test[1])

sequences = [Sequence(seq) for seq in seq_parse("data/" + FULL_SEQUENCES_FILE + ".fasta")]

for row in test:
    row[20] = -5
backtrack_matrix = BacktrackMatrixSWMSA(sequences[0], test)
# print(backtrack_matrix.s)

# print(backtrack_matrix.s)

consensus = "".join(list(map(lambda x: BASES[np.argmax(x)], pssm)))
align_smith_waterman(sequences[0], backtrack_matrix.s, consensus)
