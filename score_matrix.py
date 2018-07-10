import numpy as np


class ScoreMatrix:
    def __init__(self, seq1, seq2, sub_m, g, e, inf):
        """Generates a traceback matrix with 2 sequences and a substitution matrix.

        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.sub_m = sub_m

        self.s = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=int)
        self.v = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=int)
        self.w = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=int)

        self.init_matrices(g, e, inf)

    def init_matrices(self, g, e, inf):
        for i in range(1, len(self.s)):
            self.s[i][0] = g + (i - 1) * e
            self.v[i][0] = - inf
        for j in range(1, len(self.s[0])):
            self.s[0][j] = g + (j - 1) * e
            self.w[0][j] = - inf

        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                self.v[i][j] = max(
                    g + self.s[i][j - 1],
                    e + self.v[i][j - 1]
                )

                self.w[i][j] = max(
                    g + self.s[i - 1][j],
                    e + self.w[i - 1][j]
                )

                self.s[i][j] = max(
                    self.s[i - 1][j - 1] + self.sub_m.get_score_by_index(j - 1, self.seq1, i - 1, self.seq2),
                    self.v[i][j],
                    self.w[i][j]
                )


class ScoreMatrixSW(ScoreMatrix):
    def init_matrices(self, g, e, inf, path=[]):
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

    def recalibrate(self, path, g, e, inf):
        for index_tuple in path:
            self.s[index_tuple[0], index_tuple[1]] = 0
            self.v[index_tuple[0], index_tuple[1]] = 0
            self.w[index_tuple[0], index_tuple[1]] = 0

        self.init_matrices(g, e, inf, path)
