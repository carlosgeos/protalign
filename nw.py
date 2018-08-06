from textwrap import wrap
import re

GAP_PENALTY = -4
E_GAP_PENALTY = -1              # extended


def align_needleman_wunsch(k, sub_m, backtrack_matrix, seq1, seq2):
    s = backtrack_matrix

    i = len(seq2)
    j = len(seq1)
    align1 = ""
    align2 = ""

    score = s[i][j][0][0]
    while i > 0 or j > 0:
        direction = s[i][j][0][1]
        if i >= 0 and j >= 0 and direction == 'UL':
            align1 += seq1[j - 1]
            align2 += seq2[i - 1]
            i = i - 1
            j = j - 1
        elif direction == 'U' or (j > 0 and i == 0):
            align1 += seq1[j- 1]
            align2 += "-"
            j = j - 1
        elif direction == 'L' or (i > 0 and j == 0):
            align1 += "-"
            align2 += seq2[i - 1]
            i = i - 1


    return align1[::-1], align2[::-1], score # return flipped alignments and score
