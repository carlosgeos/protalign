import re
from functools import reduce


def dssp_parse(ifile, sequence_id="A"):
    # search for the pattern that contains the third, fourth and fifth
    # columns
    with open(ifile, 'r') as dssp:
        aa_struct_tuples = re.findall('\s+\d+\s+\d+\s' + sequence_id + '\s(.)\s\s(.)',
                                      dssp.read())
    # filter B, X and Z
    aa_struct_tuples = list(filter(lambda x: x[0] not in ['B', 'X', 'Z'],
                                   aa_struct_tuples))
    # lower case letters -> C
    aa_struct_tuples = list(map(lambda x: clean_lower_case(x),
                                aa_struct_tuples))

    return aa_struct_tuples

def clean_lower_case(tupl):
    if tupl[0].islower():
        tupl = ('C', tupl[1])
    return tupl

def create_fasta(accession, file_name, aa_struct_tuples):
    """Generates a FASTA formatted file with an accession, a file_name and
    the tuple containing the residues

    These tuples are previously extracted from the DSSP output

    """
    sequence = ''.join(map(lambda x: x[0], aa_struct_tuples))
    with open(file_name, 'a') as fasta_file:
        fasta_file.write(">" + accession + "\n")
        fasta_file.write(sequence + "\n")
