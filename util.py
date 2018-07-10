import numpy as np
import re


def seq_parse(ifile):
    """FASTA file parser. This helper function parses a FASTA formatted
    file using regular expressions

    :param ifile: input file
    :returns: all the sequences in the file
    :rtype: list

    It works as follows:

    The file is read and newline characters stripped off. A sequence
    can be discriminated from the others knowing the fact that it
    starts with the ">" separator and it ends with a string (of at
    least length 20) of upper-case letters. The lazy quantifier "*" is
    important here. We find all matches and make a first (raw) list.

    A list of tuples is then made in which each one of them contains
    the description in the first position and the amino_acid_chain in
    the second.

    Assumptions:

    - No line/sequence starts with ";"
    - Sequences do not end in "*"
    - No word > 20 chars in the description.

    FASTA file format: https://zhanglab.ccmb.med.umich.edu/FASTA/

    """
    with open(ifile, 'r') as fasta:
        data = fasta.read().replace('\n', '')
    split_sequences = re.findall('>.*?[A-Z]{20,}', data)

    for i in range(0, len(split_sequences)):
        amino_acid_chain = re.search('[A-Z]{20,}', split_sequences[i]).group(0)
        description = re.sub('[A-Z]{20,}', '', split_sequences[i])
        split_sequences[i] = (description, amino_acid_chain)

    return split_sequences


def sub_mat_parse(ifile):
    """This function parses a substitution matrix from a .bla formatted file.

    :param ifile: input file
    :returns: a matrix
    :rtype: list

    All the data is read from the file and the lines starting with '#'
    are removed. Aminoacids letters and '*' too. We are only left with
    numbers. We split this big string using the whitespace that
    separates them, thus creating a big list.

    We then pass this list to the numpy array constructor, specifying
    the integer data type and reshaping it to a 24x24 matrix.

    """

    SUB_MAT_SIZE = 24
    with open(ifile) as sub_mat:
        data = re.sub('#.*\n|\s*([A-Z]|\*)', '', sub_mat.read()).rstrip().lstrip()

        sub_mat = np.array(re.split('\s+', data), dtype=int).reshape(SUB_MAT_SIZE,SUB_MAT_SIZE)

    return sub_mat


def generate_dots(align1, align2, sub_m):
    dots = ""
    for i in range(len(align1)):
        if align1[i] == align2[i]:
            dots += ":"
        elif align1[i] == "-" or align2[i] == "-":
            dots += " "
        elif sub_m.get_score_by_letter(align1[i], align2[i]) >= 0:
            dots += "."
        else:
            dots += " "
        if len(align1) >= 80:
            align1

    return dots

def generate_lalign_output(align1, align2, dots):
    out = ""

    # split strings as a list, wrap each element to 50 characters
    a1 = re.findall('.{1,50}', align1)
    a2 = re.findall('.{1,50}', align2)
    d = re.findall('.{1,50}', dots)

    for i in range(0, len(a1)):
        # combine strings
        out += a1[i] + "\n" + d[i] + "\n" + a2[i] + "\n\n"

    return out