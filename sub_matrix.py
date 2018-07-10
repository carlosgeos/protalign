import re


class SubstitutionMatrix:
    bases = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
             'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'] # Aminoacids

    def __init__(self, substitution_mat):
        self.mat = substitution_mat

    def position_in_array(self, letter):
        return SubstitutionMatrix.bases.index(letter)

    def get_score_by_index(self, index_seq1, seq1, index_seq2, seq2):
        """Get score of letter from two indices of a sequence. To do this,
        both the index AND the sequence need to be given as argument. This
        class is not coupled to a certain pair of sequences.

        letter1_index #

        """
        letter1_index = self.position_in_array(seq1[index_seq1])
        letter2_index = self.position_in_array(seq2[index_seq2])

        return self.mat[letter1_index][letter2_index]

    def get_score_by_letter(self, letter1, letter2):
        """Get score of 2 letters.
        """
        try:
            score = self.mat[self.position_in_array(letter1)][self.position_in_array(letter2)]
        except ValueError:
            print("The letter is not an aminoacid !")
            raise
        return score

    def __str__(self):
        """String representation of the matrix in order to print it to the
        standard output.

        """
        return self.mat.__str__()
