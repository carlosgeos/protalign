from textwrap import fill


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

    def __getitem__(self, position):
        """Hook method called when using the "[]" operator to get a value at a
        certain position in the string

        :param position: string position
        :returns: the residue letter and its conformation, in a
        named dictionnary
        :rtype: string

        """
        return {"aa": self.amino_acid_chain[position],
                "conf": self.conformations[position]}

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
