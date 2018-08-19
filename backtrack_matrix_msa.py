import numpy as np
import operator

from glob_opts import GAP_PENALTY, E_GAP_PENALTY, INFINITY, BASES
g = GAP_PENALTY
e = E_GAP_PENALTY


class BacktrackMatrixSWMSA:
    def __init__(self, seq, pssm):
        """Generates a traceback matrix with 2 sequences and a substitution matrix.

        """
        self.seq = seq
        self.pssm = pssm

        # NumPy: Defining structured arrays:
        # https://docs.scipy.org/doc/numpy-1.13.0/user/basics.rec.html#structured-arrays
        self.s = np.empty((len(seq) + 1, len(pssm) + 1), dtype=object)
        self.s[0][0] = [(0, 'STOP')]
        self.init_matrices()
        self.fill_matrices()

        self.path_taken = set()

    def init_matrices(self):
        for i in range(1, len(self.s)):
            self.s[i][0] = [(0, 'STOP')]
        for j in range(1, len(self.s[0])):
            self.s[0][j] = [(0, 'STOP')]

    def fill_matrices(self, starting_i=1, starting_j=1):
        for i in range(1, len(self.s)):
            for j in range(1, len(self.s[0])):
                current_residue = self.seq[i - 1]
                gap = 20    # 21st column in PSSM
                pos = BASES.index(current_residue)
                diagonal = self.s[i - 1][j - 1][0][0] + self.pssm[j - 1][pos]
                top = self.s[i - 1][j][0][0] + self.pssm[j - 1][gap]
                left = self.s[i][j - 1][0][0] + self.pssm[j - 2][gap]

                max_option = max(diagonal, top, left, 0)

                if self.s[i][j] is None or self.s[i][j][0][1] is not "TAKEN":
                    self.s[i][j] = []
                    if(0 == max_option):
                        self.s[i][j].append((0, 'STOP', False))
                    if(diagonal == max_option):
                        self.s[i][j].append((max_option, 'UL'))
                    if(top == max_option):
                        self.s[i][j].append((max_option, 'U'))
                    if(left == max_option):
                        self.s[i][j].append((max_option, 'L'))

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
