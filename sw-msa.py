from textwrap import wrap
import re


BASES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

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


def alignSmithWaterman(pssm_m, s, seq1, g):
    path = []
    max_cell = find_max(s)

    i = max_cell['i']
    j = max_cell['j']
    final_score = 0
    align = ""
    # opening = True              # Are we opening a gap?

    while i > 0 or j > 0:
        score = s[i][j]
        score_diag = s[i - 1][j - 1]
        match_score = pssm_m[BASES.index(seq1[i - 1])][j - 1]

        if i >= 0 and j >= 0 and score == score_diag + match_score:
            path.append((i, j))
            final_score += match_score
            align += seq1[i - 1]
            i = i - 1
            j = j - 1
            opening = True      # reset flag
        elif score == s[i][j - 1] + g:
            path.append((i, j))
            final_score += g
            align += seq1[i- 1]
            j = j - 1
            opening = False     # consider next gap as extending
        elif score == s[i - 1][j] + g:
            path.append((i, j))
            final_score += g
            align += "-"
            i = i - 1
            opening = False     # consider next gap as extending

        if score == 0: break


    return align[::-1], final_score, path
