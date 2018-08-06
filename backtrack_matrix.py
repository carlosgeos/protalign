import numpy as np

GAP_PENALTY = -4
E_GAP_PENALTY = -1              # extended
INFINITY = 100000
g = GAP_PENALTY
e = E_GAP_PENALTY

class BacktrackMatrix:
    def __init__(self, seq1, seq2, sub_m):
        """Generates a traceback matrix with 2 sequences and a substitution matrix.

        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.sub_m = sub_m

        self.s = np.empty((len(seq2) + 1, len(seq1) + 1), dtype=object)
        self.v = np.empty((len(seq2) + 1, len(seq1) + 1), dtype=object)
        self.w = np.empty((len(seq2) + 1, len(seq1) + 1), dtype=object)
        self.s[0][0] = [(0,)]
        self.init_matrices(g, e)

    def init_matrices(self, g, e):
        for i in range(1, len(self.s)):
            self.s[i][0] = [(g + (i - 1) * e, 'U')]
            self.v[i][0] = - INFINITY
        for j in range(1, len(self.s[0])):
            self.s[0][j] = [(g + (j - 1) * e, 'L')]
            self.w[0][j] = - INFINITY

        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                self.v[i][j] = max(
                    g + self.s[i][j - 1][0][0],
                    e + self.v[i][j - 1]
                )

                self.w[i][j] = max(
                    g + self.s[i - 1][j][0][0],
                    e + self.w[i - 1][j]
                )
                diagonal = self.s[i - 1][j - 1][0][0] + \
                    self.sub_m.get_score_by_index(j - 1, self.seq1, i - 1, self.seq2)
                top = self.v[i][j]
                left = self.w[i][j]
                max_option = max(diagonal, top, left)

                self.s[i][j] = []

                if(diagonal == max_option):
                    self.s[i][j].append((diagonal, 'UL'))
                elif(top == max_option):
                    self.s[i][j].append((top, 'U'))
                elif(left == max_option):
                    self.s[i][j].append((left, 'L'))


class BacktrackMatrixSW(BacktrackMatrix):
    def init_matrices(self, g, e, path=[]):
        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                if (i, j) not in path:
                    self.v[i][j] = max(
                        g + self.s[i][j - 1],
                        e + self.v[i][j - 1],
                        0
                    )

                    self.w[i][j] = max(
                        g + self.s[i - 1][j],
                        e + self.w[i - 1][j],
                        0
                    )

                    self.s[i][j] = max(
                        self.s[i - 1][j - 1] + self.sub_m.get_score_by_index(j - 1, self.seq1, i - 1, self.seq2),
                        self.v[i][j],
                        self.w[i][j],
                        0
                    )

    def recalibrate(self, path, g, e):
        for index_tuple in path:
            self.s[index_tuple[0], index_tuple[1]] = 0
            self.v[index_tuple[0], index_tuple[1]] = 0
            self.w[index_tuple[0], index_tuple[1]] = 0

        self.init_matrices(g, e, path)
