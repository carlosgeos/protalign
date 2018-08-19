from score import Score
from sequence import Sequence
from backtrack_matrix import BacktrackMatrix, BacktrackMatrixSW
from nw import align_needleman_wunsch
from sw import align_smith_waterman
from util import seq_parse, sub_mat_parse
from glob_opts import WW_SEQUENCES_FILE, SUB_MATRIX, FULL_SEQUENCES_FILE


# def main():
sub_m = Score(sub_mat_parse("data/" + SUB_MATRIX + ".txt"))
ww_seqs = [Sequence(seq) for seq in seq_parse("data/" + WW_SEQUENCES_FILE + ".fasta")]

ww_seq1 = ww_seqs[0]
ww_seq2 = ww_seqs[1]
ww_seq3 = ww_seqs[2]
ww_seq4 = ww_seqs[3]
ww_seq5 = ww_seqs[4]
ww_seq6 = ww_seqs[5]

# test1 = Sequence.fromstring("ACGT")
# test2= Sequence.fromstring("ACGGCT")

# test3 = Sequence.fromstring("SLKMF")
# test4 = Sequence.fromstring("GKLKMF")

print("\nglobal aligning...\n")

k = 2                       # max number of optimal alignments
backtrack_matrix = BacktrackMatrix(ww_seq1, ww_seq2, sub_m)
align_needleman_wunsch(k, sub_m, backtrack_matrix.s, ww_seq1, ww_seq2)
backtrack_matrix = BacktrackMatrix(ww_seq3, ww_seq4, sub_m)
align_needleman_wunsch(k, sub_m, backtrack_matrix.s, ww_seq3, ww_seq4)
backtrack_matrix = BacktrackMatrix(ww_seq5, ww_seq6, sub_m)
align_needleman_wunsch(k, sub_m, backtrack_matrix.s, ww_seq5, ww_seq6)

print("\nlocal aligning...\n")

full_seqs = [Sequence(seq) for seq in seq_parse("data/" + FULL_SEQUENCES_FILE + ".fasta")]

prot1 = full_seqs[0]
prot2 = full_seqs[1]
prot3 = full_seqs[2]
prot4 = full_seqs[3]
backtrack_matrix_sw = BacktrackMatrixSW(prot1, prot2, sub_m)

l = 3
for i in range(l):
    path_taken = align_smith_waterman(sub_m, backtrack_matrix_sw.s, prot1, prot2)
    if i != l - 1:          # do not recalibrate last time
        backtrack_matrix_sw.recalibrate(path_taken)

print("\n\n...done")

# if __name__ == '__main__':
#     main()
