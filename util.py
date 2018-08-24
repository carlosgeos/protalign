import os
import contextlib
import numpy as np
import re
from glob_opts import *


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


def msa_parse(ifile):
    """Multiple alignment file parser. This helper function parses a FASTA
    formatted file using regular expressions

    :param ifile: input file
    :returns: all the alignments in the file
    :rtype: list

    It works as follows:

    The file is read and newline characters stripped off. A sequence
    can be discriminated from the others knowing the fact that it
    starts with the ">" separator and it ends with a string (of at
    least length 20) of upper-case letters and gaps ('-'). The lazy
    quantifier "*" is important here. We find all matches and make a
    first (raw) list.

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
    split_sequences = re.findall('>.*?[A-Z|-]{20,}', data)

    for i in range(0, len(split_sequences)):
        amino_acid_chain = re.search('[A-Z|-]{20,}', split_sequences[i]).group(0)
        description = re.sub('[A-Z|-]{20,}', '', split_sequences[i])
        split_sequences[i] = (description, amino_acid_chain)

    return split_sequences


def generate_dots(align1, align2, sub_m):
    """LALIGN-like display dots, the count of colons and semicolons is
    used to determine the identity and similarity percentage (just
    like in LALIGN)

    """
    dots = ""
    colons = 0
    semicolons = 0
    for i in range(len(align1)):
        if align1[i] == align2[i]:
            dots += ":"
            colons += 1
        elif align1[i] == "-" or align2[i] == "-":
            dots += " "
        elif sub_m.get_score_by_letter(align1[i], align2[i]) >= 0:
            dots += "."
            semicolons += 1
        else:
            dots += " "
        if len(align1) >= 80:
            align1

    return dots, colons, semicolons

def print_lalign_output(align1, align2, sub_m, score, seq1, seq2, align_type="NW"):
    """Mimic LALIGN's output
    """
    (dots, colons, semicolons) = generate_dots(align1, align2, sub_m)
    out = ""

    # split strings as a list, wrap each element to 60 characters
    a1 = re.findall('.{1,60}', align1)
    a2 = re.findall('.{1,60}', align2)
    d = re.findall('.{1,60}', dots)

    for i in range(0, len(a1)):
        # combine string + dots + string
        out += a1[i] + "\n" + d[i] + "\n" + a2[i] + "\n\n"

    if align_type is "NW":
        print("--- Global align ---")
    else:
        print("--- Local align ---")
    print("Opening gap penalty:\t", GAP_PENALTY)
    print("Extending gap penalty:\t", E_GAP_PENALTY)
    print(align_type, "score:\t\t", score)
    identity = "{:%}".format(colons / len(align1))
    similarity = "{:%}".format((semicolons + colons) / len(align1))
    print("Identity:\t\t", identity)
    print("Similarity:\t\t", similarity, "\n")
    print(out)


def remove_fasta_files():
    """Helper function to clean directory when we want to produce new
    FASTA files.

    """
    with contextlib.suppress(FileNotFoundError):
        os.remove(FASTA_TESTING)
    with contextlib.suppress(FileNotFoundError):
        os.remove(FASTA_TRAINING)


def parse_cath(cath_file):
    """Given a CATH_info.txt file as input, it compiles a list of (file,
    seq_ids)

    """
    list_of_f = []
    with open(cath_file, 'r') as f:
        for seq in f:
            filename = seq[:4]
            seq_id = seq[4]
            list_of_f.append((filename, seq_id))

    return list_of_f


def mapping(residue_struct_tuple):
    """Given a (residue, struct) tuple as input, modify its structure as
    follows:

    H, G, I     -> H
    E, B        -> E
    T, C, S, \s -> C

    """
    if residue_struct_tuple[1] in ['H', 'G', 'I']:
        residue_struct_tuple = (residue_struct_tuple[0], 'H')
    elif residue_struct_tuple[1] in ['E', 'B']:
        residue_struct_tuple = (residue_struct_tuple[0], 'E')
    else:
        residue_struct_tuple = (residue_struct_tuple[0], 'C')

    return residue_struct_tuple


def clean_lower_case(tupl):
    """From the documentation:

    If this (AA) letter is lower case this means this is a cysteine
    that form a sulfur bridge with the other amino acid in this column
    with the same lower case letter.

    We simply set it to C.

    """
    if tupl[0].islower():
        tupl = ('C', tupl[1])
    return tupl


def dssp_parse(file_id, sequence_id="A"):
    """Parses a DSSP formatted file and returns a list of the form:

    ("ARNDCQEGHILKMFPSTWYV", "CEHCEHCEHCEHCEHCEHCEH")

    """
    if CATH_FILE == TRAINING_CATH_FILE:
        full_path = DSSP_PATH + file_id + ".dssp"
    else:
        full_path = DSSP_TEST_PATH + file_id + ".dssp"

    # search for the pattern that contains the third, fourth and fifth
    # columns
    with open(full_path, 'r') as dssp:
        aa_struct_tuples = re.findall('\s+\d+\s+\d+\s' + sequence_id + '\s(.)\s\s(.)',
                                      dssp.read())
    # filter B, X and Z
    aa_struct_tuples = list(filter(lambda x: x[0] not in ['B', 'X', 'Z'],
                                   aa_struct_tuples))
    # lower case letters -> C
    aa_struct_tuples = list(map(lambda x: clean_lower_case(x),
                                aa_struct_tuples))

    # apply mapping
    aa_struct_tuples = list(map(lambda x: mapping(x), aa_struct_tuples))

    sequence = ''.join(map(lambda x: x[0], aa_struct_tuples))
    structures = ''.join(map(lambda x: x[1], aa_struct_tuples))

    return (sequence, structures)


def create_fasta(accession, sequence, structure):
    """Generates a FASTA formatted file with an accession, a file_name and
    the tuple containing the residues

    These tuples are previously extracted from the DSSP output

    The only thing that is different is that its structure is also
    included in the next line

    """
    if CATH_FILE == TRAINING_CATH_FILE:
        file_name = FASTA_TRAINING
    else:
        file_name = FASTA_TESTING

    with open(file_name, 'a') as fasta_file:
        fasta_file.write(">" + accession + "\n")
        fasta_file.write(sequence + "\n")
        fasta_file.write(structure + "\n")
