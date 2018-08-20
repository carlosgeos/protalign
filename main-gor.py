from util import remove_fasta_files, parse_cath, dssp_parse, create_fasta
from glob_opts import *


# def main():
remove_fasta_files()

list_of_files = parse_cath(CATH_FILE)
for input_dssp in list_of_files:
    if CATH_FILE == TRAINING_CATH_FILE:
        aa_struct_tuples = dssp_parse(DSSP_PATH + input_dssp[0] + ".dssp", input_dssp[1])
        create_fasta(input_dssp[0] + input_dssp[1], FASTA_TRAINING, aa_struct_tuples)
    else:
        aa_struct_tuples = dssp_parse(DSSP_TEST_PATH + input_dssp[0] + ".dssp", input_dssp[1])
        create_fasta(input_dssp[0] + input_dssp[1], FASTA_TESTING, aa_struct_tuples)

# From this poin on we've got our FASTA file



# if __name__ == '__main__':
#     main()
