from textwrap import wrap
import re

def find_max(s):
    max_cell = {
        'i': 0,
        'j': 0,
        'max_value': 0
    }
    for i in range(len(s)):
        for j in range(len(s[0])):
            current_value = s[i][j]
            if current_value > max_cell['max_value']:
                max_cell['i'] = i
                max_cell['j'] = j
                max_cell['max_value'] = current_value

    return max_cell


def align_smith_waterman(sub_m, svw, seq1, seq2, g, e):
    s = svw[0]; v = svw[1]; w = svw[2] # unpacking
    path = []
    max_cell = find_max(s)

    i = max_cell['i']
    j = max_cell['j']
    final_score = 0
    align1 = ""
    align2 = ""
    opening = True              # Are we opening a gap?

    while i > 0 or j > 0:
        score = s[i][j]
        score_diag = s[i - 1][j - 1]
        match_score = sub_m.get_score_by_index(j - 1, seq1, i - 1, seq2)

        if i >= 0 and j >= 0 and score == score_diag + match_score:
            path.append((i, j))
            final_score += match_score
            align1 += seq1[j - 1]
            align2 += seq2[i - 1]
            i = i - 1
            j = j - 1
            opening = True      # reset flag

        elif score == v[i][j]:
            path.append((i, j))
            final_score += g if opening else e
            align1 += seq1[j- 1]
            align2 += "-"
            j = j - 1
            opening = False     # consider next gap as extending
        elif score == w[i][j]:
            path.append((i, j))
            final_score += g if opening else e
            align1 += "-"
            align2 += seq2[i - 1]
            i = i - 1
            opening = False     # consider next gap as extending

        if score == 0: break


    return align1[::-1], align2[::-1], final_score, path
