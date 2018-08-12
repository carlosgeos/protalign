from util import print_lalign_output


def align_needleman_wunsch(sub_m, backtrack_matrix, seq1, seq2,
                           i=None, j=None, align1="", align2="",
                           extending=False):
    s = backtrack_matrix

    if i is None:
        i = len(seq2)
    if j is None:
        j = len(seq1)

    # score is located in the last cell, first tuple, first item
    score = s[len(seq2)][len(seq1)][0][0]

    while (i > 0 or j > 0):
        direction = s[i][j][0][1]

        # means extending
        if s[i][j][0][2]:
            # remove the UL option in next cell in path, since we are
            # extending and either L or U should follow, but never UL
            if direction == 'L':
                s[i][j - 1] = list(filter(lambda x: x[1] != 'UL', s[i][j - 1]))
            elif direction == 'U':
                s[i - 1][j] = list(filter(lambda x: x[1] != 'UL', s[i - 1][j]))

        # check for more path possiblities (optimal)
        if len(s[i][j]) > 1:
            # pop path arrow we are exploring below and call function
            # the next one
            s[i][j].pop(0)
            align_needleman_wunsch(sub_m, s, seq1, seq2, i, j, align1, align2)

        if i >= 0 and j >= 0 and direction == 'UL' and not extending:
            align1 += seq1[j - 1]
            align2 += seq2[i - 1]
            i = i - 1
            j = j - 1
        elif direction == 'L':
            align1 += seq1[j - 1]
            align2 += "-"
            j = j - 1
        elif direction == 'U':
            align1 += "-"
            align2 += seq2[i - 1]
            i = i - 1

    print_lalign_output(align1[::-1], align2[::-1], sub_m, score, seq1, seq2)
