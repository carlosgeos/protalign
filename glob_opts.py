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
