from sub_matrix import SubstitutionMatrix
from sequence import Sequence
from score_matrix import ScoreMatrix, ScoreMatrixSW
from nw import alignNeedlemanWunsch
from sw import alignSmithWaterman

from util import seq_parse, sub_mat_parse, generate_dots, generate_lalign_output


SEQUENCE_FILE = "protein-sequences"
SUB_MATRIX = "blosum62"
GAP_PENALTY = -12
EXTENDED_GAP_PENALTY = -2
INFINITY = -1000

def main():
    seq_all = seq_parse("data/" + SEQUENCE_FILE + ".fasta")
    sequences = [Sequence(seq) for seq in seq_all]

    seq1 = sequences[0]; seq2 = sequences[1]

    #seq1 = Sequence.fromstring("MGGETFA")
    #seq2 = Sequence.fromstring("GGVTTF")

    sub_m = SubstitutionMatrix(sub_mat_parse("data/" + SUB_MATRIX + ".txt"))

    score_matrix = ScoreMatrix(seq1, seq2, sub_m, GAP_PENALTY,
                               EXTENDED_GAP_PENALTY, INFINITY)
    svw = (score_matrix.s, score_matrix.v, score_matrix.w) # Three packed matrices

    k = 1
    align1, align2, final_score = alignNeedlemanWunsch(k, sub_m, svw, seq1, seq2,
                                                       GAP_PENALTY, EXTENDED_GAP_PENALTY)
    dots = generate_dots(align1, align2, sub_m)
    ll_output = generate_lalign_output(align1, align2, dots)


    print("NW align")
    print("Similarity score:", final_score)
    print(ll_output)

    print("---------------------------------------------------")
    print("SW align")
    # SW
    score_matrix_sw = ScoreMatrixSW(seq1, seq2, sub_m, GAP_PENALTY,
                                    EXTENDED_GAP_PENALTY, INFINITY)

    svw = (score_matrix_sw.s, score_matrix_sw.v, score_matrix_sw.w)

    l = 3
    for i in range(0, l):
        align1, align2, final_score, path = alignSmithWaterman(sub_m, svw, seq1, seq2,
                                                               GAP_PENALTY, EXTENDED_GAP_PENALTY)
        score_matrix_sw.recalibrate(path, GAP_PENALTY, EXTENDED_GAP_PENALTY, INFINITY)

        print("Local alignment number:", i + 1)
        dots = generate_dots(align1, align2, sub_m)
        ll_output = generate_lalign_output(align1, align2, dots)
        print("Similarity score:", final_score)
        print(ll_output)


if __name__ == '__main__':
    main()
