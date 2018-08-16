import numpy as np
import operator

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
        self.w = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=struct_array)
        self.s[0][0] = [(0, 'STOP', False)]
        self.init_matrices()
        self.fill_matrices()

    def init_matrices(self):
        for i in range(1, len(self.s)):
            self.s[i][0] = [(g + (i - 1) * e, 'U', False)]
            self.v[i][0] = (- INFINITY, False)
        for j in range(1, len(self.s[0])):
            self.s[0][j] = [(g + (j - 1) * e, 'L', False)]
            self.w[0][j] = (- INFINITY, False)

    def fill_matrices(self):
        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                max_left = max(g + self.s[i][j - 1][0][0],
                               e + self.v[i][j - 1]['score'])

                self.v[i][j]['score'] = max_left

                if e + self.v[i][j - 1]['score'] == max_left:
                    self.v[i][j]['extending?'] = True
                else:
                    self.v[i][j]['extending?'] = False

                max_top = max(g + self.s[i - 1][j][0][0],
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
                    self.s[i][j].append((max_option, 'UL', False))
                if(top == max_option):
                    self.s[i][j].append((max_option, 'U', self.w[i][j]['extending?']))
                if(left == max_option):
                    self.s[i][j].append((max_option, 'L', self.v[i][j]['extending?']))


class BacktrackMatrixSW(BacktrackMatrix):
    def __init__(self, seq1, seq2, sub_m):
        """Generates a traceback matrix with 2 sequences and a substitution matrix.

        """
        super().__init__(seq1, seq2, sub_m)
        self.path_taken = set()

    def init_matrices(self):
        for i in range(1, len(self.s)):
            self.s[i][0] = [(0, 'STOP', False)]
        for j in range(1, len(self.s[0])):
            self.s[0][j] = [(0, 'STOP', False)]

    def fill_matrices(self, starting_i=1, starting_j=1):
        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                max_left = max(g + self.s[i][j - 1][0][0],
                               e + self.v[i][j - 1]['score'],
                               0)

                self.v[i][j]['score'] = max_left

                if e + self.v[i][j - 1]['score'] == max_left:
                    self.v[i][j]['extending?'] = True
                else:
                    self.v[i][j]['extending?'] = False

                max_top = max(g + self.s[i - 1][j][0][0],
                              e + self.w[i - 1][j]['score'],
                              0)

                self.w[i][j]['score'] = max_top

                if e + self.w[i - 1][j]['score'] == max_top:
                    self.w[i][j]['extending?'] = True
                else:
                    self.w[i][j]['extending?'] = False

                diagonal = self.s[i - 1][j - 1][0][0] + \
                    self.sub_m.get_score_by_index(j - 1, self.seq1, i - 1, self.seq2)
                left = self.v[i][j]['score']
                top = self.w[i][j]['score']
                max_option = max(diagonal, top, left, 0)

                if self.s[i][j] is None or self.s[i][j][0][1] is not "TAKEN":
                    self.s[i][j] = []
                    if(0 == max_option):
                        self.s[i][j].append((0, 'STOP', False))
                    if(diagonal == max_option):
                        self.s[i][j].append((max_option, 'UL', False))
                    if(top == max_option):
                        self.s[i][j].append((max_option, 'U', self.w[i][j]['extending?']))
                    if(left == max_option):
                        self.s[i][j].append((max_option, 'L', self.v[i][j]['extending?']))

    def recalibrate(self, path_taken):
        self.path_taken = self.path_taken.union(path_taken)

        for cell_tuple in self.path_taken:
            self.s[cell_tuple[0], cell_tuple[1]] = [(0, 'TAKEN', False)]

        # we sort the path, by rows and then tuples, to obtain the
        # minimum i and j with which we need to start recalibrating
        # the matrix
        sorted_path = sorted(list(self.path_taken), key=operator.itemgetter(0, 1))
        starting_i = sorted_path[0][0]
        starting_j = sorted_path[0][1]
        self.init_matrices()
        self.fill_matrices(starting_i, starting_j)
