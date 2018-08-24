from util import remove_fasta_files, parse_cath, dssp_parse, create_fasta
from sequence_struct import SequenceWithStructure
from glob_opts import *


# def main():
remove_fasta_files()

sequences = []
list_of_files = parse_cath(CATH_FILE)
for input_dssp in list_of_files:
    accession = input_dssp[0] + input_dssp[1]
    sequence, structure = dssp_parse(input_dssp[0], input_dssp[1])
    sequences.append(SequenceWithStructure([accession, sequence, structure]))
    create_fasta(input_dssp[0] + input_dssp[1], sequence, structure)

# From this point on we've got a FASTA file written and the sequences
# in memory




# if __name__ == '__main__':
#     main()
