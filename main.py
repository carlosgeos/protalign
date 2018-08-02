from score import Score
from sequence import Sequence
from backtrack_matrix import BacktrackMatrix, BacktrackMatrixSW
from nw import align_needleman_wunsch
from sw import align_smith_waterman

from util import seq_parse, sub_mat_parse, generate_dots, generate_lalign_output


SEQUENCE_FILE = "WW-sequence"
SUB_MATRIX = "blosum50"
GAP_PENALTY = -4
E_GAP_PENALTY = -1              # extended
INFINITY = -100000

def main():
    seq_all = seq_parse("data/" + SEQUENCE_FILE + ".fasta")
    sequences = [Sequence(seq) for seq in seq_all]

    seq1 = sequences[1]; seq2 = sequences[2]

    #seq1 = Sequence.fromstring("THISLINE")
    #seq2 = Sequence.fromstring("ISALIGNED")

    sub_m = Score(sub_mat_parse("data/" + SUB_MATRIX + ".txt"))

    score_matrix = BacktrackMatrix(seq1, seq2, sub_m, GAP_PENALTY,
                                   E_GAP_PENALTY, INFINITY)
    svw = score_matrix.packed() # Three packed matrices in a tuple

    k = 1
    align1, align2, final_score = align_needleman_wunsch(k, sub_m, svw, seq1, seq2,
                                                         GAP_PENALTY, E_GAP_PENALTY)
    dots = generate_dots(align1, align2, sub_m)
    ll_output = generate_lalign_output(align1, align2, dots)


    print("--- Global align ---")
    print("Opening gap penalty:\t", GAP_PENALTY)
    print("Extending gap penalty:\t", E_GAP_PENALTY)
    print("N-W score:\t\t", final_score)
    identity = "{:%}".format(ll_output["colons"] / max(len(seq1), len(seq2)))
    similarity = "{:%}".format((ll_output["semicolons"] + ll_output["colons"])
                               / max(len(seq1), len(seq2)))
    print("Identity:\t\t", identity)
    print("Similarity:\t\t", similarity, "\n")
    print(ll_output["aligned_sequences"])

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
