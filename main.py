from score import Score
from sequence import Sequence
from backtrack_matrix import BacktrackMatrix, BacktrackMatrixSW
from nw import align_needleman_wunsch
from sw import align_smith_waterman
from util import seq_parse, sub_mat_parse
from glob_opts import SEQUENCE_FILE, SUB_MATRIX


def main():
    seq_all = seq_parse("data/" + SEQUENCE_FILE + ".fasta")
    sequences = [Sequence(seq) for seq in seq_all]

    #seq2 = sequences[2]; seq1 = sequences[3]
    # print(seq1)
    # print(seq2)

    seq1 = Sequence.fromstring("ACGT")
    seq2 = Sequence.fromstring("ACGGCT")

    # seq1 = Sequence.fromstring("SLKMF")
    # seq2 = Sequence.fromstring("GKLKMF")

    sub_m = Score(sub_mat_parse("data/" + SUB_MATRIX + ".txt"))

    backtrack_matrix = BacktrackMatrix(seq1, seq2, sub_m)
    align_needleman_wunsch(sub_m, backtrack_matrix.s, seq1, seq2)


    # print("---------------------------------------------------")
    # print("SW align")
    # SW
    # score_matrix_sw = BacktrackMatrixSW(seq1, seq2, sub_m, GAP_PENALTY,
    #E_GAP_PENALTY, INFINITY)

    # svw = (score_matrix_sw.s, score_matrix_sw.v, score_matrix_sw.w)

    # l = 3
    # for i in range(0, l):
    #     align1, align2, final_score, path = align_smith_waterman(sub_m, svw, seq1, seq2,
    #                                                              GAP_PENALTY, E_GAP_PENALTY)
    #     score_matrix_sw.recalibrate(path, GAP_PENALTY, E_GAP_PENALTY, INFINITY)

    #     print("Local alignment number:", i + 1)
    #     dots = generate_dots(align1, align2, sub_m)
    #     ll_output = generate_lalign_output(align1, align2, dots)
    #     print("Similarity score:", final_score)
    #     print(ll_output)


if __name__ == '__main__':
    main()
