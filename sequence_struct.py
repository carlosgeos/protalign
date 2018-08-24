from textwrap import fill
from glob_opts import BASES


class SequenceWithStructure:
    def __init__(self, info):
        """Amino-acid sequence abstract data type.

        :param info: tuple containing description and aa-chain.
        :returns: a new Sequence object
        :rtype: Sequence

        The param info is supposed to be in a structured manner
        (obtained from the sequence parser). The description must be
        in the first position and the aa-chain in the second.

        """
        self.description = info[0]
        self.amino_acid_chain = info[1]
        self.conformations = info[2]
        self.prediction = ""

    @classmethod
    def fromstring(cls, string, structure):
        """This 'constructor' overload allows creating a sequence from a
        simple string. The information does not need to be in a
        FASTA-like format

        :param cls: our class
        :param string:
        :returns: a new object created with the default constructor
        :rtype: Sequence

        """

        info = ["No description", string, structure]
        return cls(info)

    def conf(self):
        """simply return the secondary structure
        """
        return self.conformations

    def neighbours(self, i):
        """Returns the R_j+m (m taking values from -8 to +8) indices of the
        neighbours of aminoacid at position j (index)

        """
        # boundaries
        left_b = max(0, i - 8)
        right_b = min(1 + i + 8, len(self.amino_acid_chain) - 1)
        all = list(self.amino_acid_chain[left_b:right_b])
        if right_b != len(self.amino_acid_chain) - 1:
            all.remove(self["aa", i])
        all = list(map(lambda aa: BASES.index(aa), all)) # return indices
        return all

    def __getitem__(self, selection):
        """Hook method called when using the "[]" operator to get a value at a
        certain position in the string

        :param selection: aa or conf and position
        :returns: the residue letter and its conformation
        :rtype: str

        """
        if selection[0] == "aa":
            return self.amino_acid_chain[selection[1]]
        elif selection[0] == "conf":
            return self.conformations[selection[1]]
        elif selection[0] == "pred":
            return self.prediction[selection[1]]

    def __len__(self):
        """Hook method. len() on objects with Sequence type will call this.

        :returns: length of the aa-chain
        :rtype: int

        """

        return len(self.amino_acid_chain)

    def __str__(self):
        """FASTA format representation of the sequence and its secondary
        structure

        It wraps the amino acid chain to 80 characters, which is the
        standard for the FASTA format. The description is left
        unwrapped.

        :returns:
        :rtype: string

        """

        aa = fill(self.amino_acid_chain + "\n", 80)
        conf = fill(self.conformations + "\n", 80)
        string = self.description + "\n" + aa + "\n" + conf

        return string
