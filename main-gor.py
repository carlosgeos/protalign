from util import remove_fasta_files, parse_cath, dssp_parse, create_fasta
from sequence_struct import SequenceWithStructure
from math import log
import numpy as np
from glob_opts import *
from collections import Counter
from functools import reduce
from gor import *


# def main():
remove_fasta_files()

sequences = []
list_of_files = parse_cath(CATH_FILE)
proteins_processed = set()      # limit to 700

print("\nParsing DSSP files...", end="")
for input_dssp in list_of_files:
    # we only use the 700 first proteins to train GOR
    accession = input_dssp[0] + input_dssp[1]
    sequence, structure = dssp_parse(input_dssp[0], input_dssp[1])
    sequences.append(SequenceWithStructure([">" + accession, sequence, structure]))
    create_fasta(input_dssp[0] + input_dssp[1], sequence, structure)
    proteins_processed = proteins_processed.union({input_dssp[0]})
    if len(proteins_processed) >= 700:
        break

print("done\n")



print("\n\nFilling counters...", end="")
counters = np.zeros((len(CONFORMATIONS), len(BASES), len(BASES) + 1), dtype=int)
counters = fill_counters(sequences, counters)
print("done\n")


# Filter training and prediction
files_to_predict = list(filter(lambda x: x[0] not in proteins_processed, list_of_files))
seqs_to_predict = []
for input_dssp in files_to_predict:
    accession = input_dssp[0] + input_dssp[1]
    sequence, structure = dssp_parse(input_dssp[0], input_dssp[1])
    seqs_to_predict.append(SequenceWithStructure([">" + accession, sequence, structure]))


print("\n\nPredicting secondary structure for", len(seqs_to_predict), "sequences...", end="")
for seq in seqs_to_predict:
    prediction = gor(seq, counters)
    seq.prediction = prediction # save it

print("done\n")

print("First ten predictions:\n")
# Print 10 sequences
for i in range(10):
    seq = seqs_to_predict[i]
    print(seq.description)
    print("DSSP:")
    print(seq.conformations, "\n")
    print("PRED:")
    print(seq.prediction)
    print("\n")

# now we compare seq.conformations and seq.prediction, both are in the Sequence object
q3s = []
for seq in seqs_to_predict:
    total_confs = 0
    total_correct_confs = 0
    for i in range(len(seq)):
        total_confs += 1
        if seq["conf", i] == seq["pred", i]:
            total_correct_confs += 1
    q3s.append(total_correct_confs / total_confs)

overall_q3 = reduce(lambda x, y: x + y, q3s) / len(q3s)

print("The overall Q3 score is:", overall_q3, "\n\n")

# if __name__ == '__main__':
#     main()
