from util import sub_mat_parse, print_lalign_output
from glob_opts import BASES, SUB_MATRIX
from score import Score
import numpy as np


def find_max(s):
    """Finds the max value in the matrix, together with its indeces

    """
    max_cell = {
        'i': 0,
        'j': 0,
        'max_value': 0
    }
    for i in range(len(s)):
        for j in range(len(s[0])):
            current_value = s[i][j][0][0]
            if current_value > max_cell['max_value']:
                max_cell['i'] = i
                max_cell['j'] = j
                max_cell['max_value'] = current_value

    return max_cell


def align_smith_waterman(seq, s, consensus):
    """Starting at the max value in the matrix, bactrack to a cell whose
    score is 0, and align the partial sequences.

    """
    # this substitution matrix is only used for the dot display.
    sub_m = Score(sub_mat_parse("data/" + SUB_MATRIX + ".txt"))

    max_cell = find_max(s)

    i = max_cell['i']
    j = max_cell['j']
    score = s[i][j][0][0]
    align1 = ""
    align2 = ""
    path = set()

    while (i > 0 or j > 0):
        print("SW !!!")
        print(i)
        print(j)
        path.add((i, j))     # this path will be set to 0 when recalibrating
        cell_score = s[i][j][0][0]
        print("cell_score: ", cell_score)
        direction = s[i][j][0][1]
        # means extending
        # if s[i][j][0][2]:
        #     # remove the UL option in next cell in path, since we are
        #     # extending and either L or U should follow, but never UL
        #     if direction == 'L':
        #         if len(s[i][j - 1]) > 1:
        #             s[i][j - 1] = list(filter(lambda x: x[1] != 'UL', s[i][j - 1]))
        #     elif direction == 'U':
        #         if len(s[i - 1][j]) > 1:
        #             s[i - 1][j] = list(filter(lambda x: x[1] != 'UL', s[i - 1][j]))

        if i >= 0 and j >= 0 and direction == 'UL':
            align1 += consensus[j - 1]
            align2 += seq[i - 1]
            i = i - 1
            j = j - 1
        elif direction == 'L':
            align1 += consensus[j - 1]
            align2 += "-"
            j = j - 1
        elif direction == 'U':
            align1 += "-"
            align2 += seq[i - 1]
            i = i - 1

        if cell_score <= 0: break

    print_lalign_output(align1[::-1], align2[::-1], sub_m, score, seq, consensus, "SW")

    return path                 # for recalibration
