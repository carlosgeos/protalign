import os
import contextlib
from parser import dssp_parse, create_fasta

DSSP_PATH = "data/dssp/"
CATH_FILE = "data/CATH_info.txt"

FASTA_TRAINING = "trainingset.fasta"
FASTA_TESTING = "testingset.fasta"

def clean_fasta_files():
    # we are appending so we want nothing there
    with contextlib.suppress(FileNotFoundError):
        os.remove(FASTA_TRAINING)
        os.remove(FASTA_TESTING)


def parse_cath(cath_file):
    list_of_f = []
    with open(cath_file, 'r') as f:
        for seq in f:
            filename = seq[:4]
            seq_id = seq[4]
            list_of_f.append((filename, seq_id))

    return list_of_f

def main():
    clean_fasta_files()

    list_of_files = parse_cath(CATH_FILE)
    # eso = dssp_parse("data/dssp_test/1AVA.dssp")
    for ifile in list_of_files:
        aa_struct_tuples = dssp_parse(DSSP_PATH + ifile[0] + ".dssp", ifile[1])
        create_fasta(ifile[0] + ifile[1], FASTA_TRAINING, aa_struct_tuples)


if __name__ == '__main__':
    main()
