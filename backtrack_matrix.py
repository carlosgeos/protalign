import numpy as np

from glob_opts import GAP_PENALTY, E_GAP_PENALTY, INFINITY
g = GAP_PENALTY
e = E_GAP_PENALTY


class BacktrackMatrix:
    def __init__(self, seq1, seq2, sub_m):
        """Generates a traceback matrix with 2 sequences and a substitution matrix.

        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.sub_m = sub_m

        # NumPy: Defining structured arrays:
        # https://docs.scipy.org/doc/numpy-1.13.0/user/basics.rec.html#structured-arrays
        struct_array = {'names': ['score', 'extending?'],
                        'formats': [np.int32, 'b']}
        self.s = np.empty((len(seq2) + 1, len(seq1) + 1), dtype=object)
        self.v = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=struct_array)
        self.w = np.empty((len(seq2) + 1, len(seq1) + 1), dtype=struct_array)
        self.s[0][0] = [(0, 'STOP', False)]
        self.init_matrices(g, e)

    def init_matrices(self, g, e):
        for i in range(1, len(self.s)):
            self.s[i][0] = [(g + (i - 1) * e, 'U', False)]
            self.v[i][0] = (- INFINITY, False)
        for j in range(1, len(self.s[0])):
            self.s[0][j] = [(g + (j - 1) * e, 'L', False)]
            self.w[0][j]['score'] = - INFINITY

        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                max_left = max(g + self.s[i][j - 1][0][0],
                               e + self.v[i][j - 1]['score'])

                self.v[i][j]['score'] = max_left

                if e + self.v[i][j - 1]['score'] == max_left:
                    self.v[i][j]['extending?'] = True
                else:
                    self.v[i][j]['extending?'] = False

                max_top = max(
                    g + self.s[i - 1][j][0][0],
                    e + self.w[i - 1][j]['score'])

                self.w[i][j]['score'] = max_top

                if e + self.w[i - 1][j]['score'] == max_top:
                    self.w[i][j]['extending?'] = True
                else:
                    self.w[i][j]['extending?'] = False

                diagonal = self.s[i - 1][j - 1][0][0] + \
                    self.sub_m.get_score_by_index(j - 1, self.seq1, i - 1, self.seq2)
                left = self.v[i][j]['score']
                top = self.w[i][j]['score']
                max_option = max(diagonal, top, left)

                self.s[i][j] = []

                if(diagonal == max_option):
                    self.s[i][j].append((diagonal, 'UL', False))
                if(top == max_option):
                    self.s[i][j].append((top, 'U', self.w[i][j]['extending?']))
                if(left == max_option):
                    self.s[i][j].append((left, 'L', False))


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
