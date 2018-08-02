from textwrap import wrap
import re

def align_needleman_wunsch(k, sub_m, svw, seq1, seq2, g, e):
    s = svw[0]; v = svw[1]; w = svw[2] # unpacking

    i = len(seq2)
    j = len(seq1)
    final_score = 0
    align1 = ""
    align2 = ""
    opening = True              # Are we opening a gap?

    while i > 0 or j > 0:
        score = s[i][j]
        score_diag = s[i - 1][j - 1]
        match_score = sub_m.get_score_by_index(j - 1, seq1, i - 1, seq2)

        if i >= 0 and j >= 0 and score == score_diag + match_score:
            final_score += match_score
            align1 += seq1[j - 1]
            align2 += seq2[i - 1]
            i = i - 1
            j = j - 1
            opening = True      # reset flag
        elif score == v[i][j] or (j > 0 and i == 0):
            final_score += g if opening else e
            align1 += seq1[j- 1]
            align2 += "-"
            j = j - 1
            opening = False     # consider next gap as extending
        elif score == w[i][j] or (i > 0 and j == 0):
            final_score += g if opening else e
            align1 += "-"
            align2 += seq2[i - 1]
            i = i - 1
            opening = False     # consider next gap as extending


    return align1[::-1], align2[::-1], final_score # return flipped alignments and score
