# Pairwise alignment
WW_SEQUENCES_FILE = "WW-sequence"
FULL_SEQUENCES_FILE = "protein-sequences"
SUB_MATRIX = "blosum50"
GAP_PENALTY = -2
E_GAP_PENALTY = 0               # extended
INFINITY = 100000


# Multiple alignment
ALIGNMENTS_FILE_NAMES = {"muscle": "data/msaresults-muscle.fasta",
                         "clustal-omega": "data/msaresults-clustalo.fasta",
                         "tcoffee": "data/msaresults-tcoffee.fasta"}

BASES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']


# GOR
BASES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
f_sr_i = len(BASES)         # index to store f_sr

CATH_FILE = "data/CATH_info.txt"  # change !

DSSP_PATH = "data/dssp/"
DSSP_TEST_PATH = "data/dssp_test/"
TEST_CATH_FILE = "data/CATH_info_test.txt"
TRAINING_CATH_FILE = "data/CATH_info.txt"

FASTA_TRAINING = "trainingset.fasta"
FASTA_TESTING = "testingset.fasta"

CONFORMATIONS = ['C', 'E', 'H']
CONFORMATIONS_INDICES = [0, 1, 2]
