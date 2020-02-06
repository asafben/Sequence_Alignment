# *** Imports *** #
import argparse
import numpy as np
from itertools import groupby
import csv

# *** Constants *** #
GAP = '-'
UP = 1
LEFT = 2
DIAGONAL = 0
GLOBAL = "global"
LOCAL = "local"
BUFF_SIZE = 50

# *** Functions *** #
def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def print_alignment(seq1, seq2, align_type, score):
    buff_repeats = (max(len(seq1), len(seq2)) // BUFF_SIZE) + 1  # Round up.
    for i in range(0, buff_repeats):
        print seq1[i * BUFF_SIZE : (i + 1) * BUFF_SIZE]
        print seq2[i * BUFF_SIZE : (i + 1) * BUFF_SIZE] + '\n'

    print align_type + ':' + str(score)


def traceback(dynamic_table, seq1, seq2, align_type, max_i, max_j):
    # Both seq after alignment.
    ali_seq1 = ""
    ali_seq2 = ""

    while max_i > 0 or max_j > 0:
        # additional stop condition, in case of local alignment.
        if dynamic_table[max_i][max_j][0] == 0 and align_type == LOCAL:
            break

        if (max_i > 0) and dynamic_table[max_i][max_j][1] == UP:
            ali_seq1 = seq1[max_i - 1] + ali_seq1
            ali_seq2 = "-" + ali_seq2
            max_i -= 1

        elif dynamic_table[max_i][max_j][1] == LEFT:
            ali_seq1 = "-" + ali_seq1
            ali_seq2 = seq2[max_j - 1] + ali_seq2
            max_j -= 1

        elif (max_i > 0) and (max_j > 0) and dynamic_table[max_i][max_j][1] == DIAGONAL:
            ali_seq1 = seq1[max_i - 1] + ali_seq1
            ali_seq2 = seq2[max_j - 1] + ali_seq2
            max_i -= 1
            max_j -= 1

    return ali_seq1, ali_seq2


# Update the rest of the table values, from the top left corner down - col by col.
def fill_inner_table(dynamic_table, base_to_base_score_map, seq1, seq2, seq1_len, seq2_len, align_type):
    for i in range(1, seq1_len):      # Row num
        for j in range(1, seq2_len):  # Col num

            # Update score.
            match_gain = dynamic_table[i-1][j-1][0] + base_to_base_score_map[seq1[i-1]][seq2[j-1]]  # Diagonal - 0
            gap_gain_seq_1 = dynamic_table[i-1][j][0] + base_to_base_score_map[seq1[i-1]][GAP]  # Up - 1
            gap_gain_seq_2 = dynamic_table[i][j-1][0] + base_to_base_score_map[GAP][seq2[j-1]]  # LEFT - 2

            gains_arr = np.array([match_gain, gap_gain_seq_1, gap_gain_seq_2])

            dynamic_table[i][j][0] = np.max(gains_arr)     # Current cell new value.
            dynamic_table[i][j][1] = np.argmax(gains_arr)  # Update previous cell direction.

            # In case of local alignment, no account to previous penalties.
            if dynamic_table[i][j][0] < 0 and align_type == LOCAL:
                dynamic_table[i][j][0] = 0

'''
    Each matrix cell (val,i,j,direction) holds:
        val - Calculated score value.
        direction - direction of the previous cell:
                    (i)   Up       - (i-1, j)   = 0
                    (ii)  LEFT     - (i, j-1)   = 1
                    (iii) Diagonal - (i-1, j-1) = 2
'''


def global_seq_align(dynamic_table, base_to_base_score_map, seq1, seq2, seq1_len, seq2_len):

    # Initialize first row with the gap penalty.
    for i in range(1, seq2_len): # Iterate over num of cols.
        dynamic_table[0][i][0] = dynamic_table[0][i-1][0] + base_to_base_score_map[GAP][seq2[i-1]]  # Update score.
        dynamic_table[0][i][1] = LEFT                                     # Update previous cell direction.

    # Initialize first col with the gap penalty.
    for i in range(1, seq1_len): # Iterate over num of rows.
        dynamic_table[i][0][0] = dynamic_table[i-1][0][0] + base_to_base_score_map[seq1[i-1]][GAP]  # Update score.
        dynamic_table[i][0][1] = UP                                     # Update previous cell direction.

    fill_inner_table(dynamic_table, base_to_base_score_map, seq1, seq2, seq1_len, seq2_len, GLOBAL)
    seq1, seq2 = traceback(dynamic_table, seq1, seq2, GLOBAL, seq1_len - 1, seq2_len - 1)

    return seq1, seq2, dynamic_table[seq1_len - 1, seq2_len - 1][0]


def local_seq_align(dynamic_table, base_to_base_score_map, seq1, seq2, seq1_len, seq2_len):
    # Initialize first row with the gap penalty.
    for i in range(1, seq2_len):  # Iterate over num of cols.
        dynamic_table[0][i][0] = 0     # Update score.
        dynamic_table[0][i][1] = LEFT  # Update previous cell direction.

    # Initialize first col with the gap penalty.
    for i in range(1, seq1_len):  # Iterate over num of rows.
        dynamic_table[i][0][0] = 0   # Update score.
        dynamic_table[i][0][1] = UP  # Update previous cell direction.

    fill_inner_table(dynamic_table, base_to_base_score_map, seq1, seq2, seq1_len, seq2_len, LOCAL)
    start_i, start_j = np.unravel_index(np.argmax(dynamic_table[:, :, 0]), dynamic_table.shape[0:2])
    seq1, seq2 = traceback(dynamic_table, seq1, seq2, LOCAL, start_i, start_j)

    return seq1, seq2, dynamic_table[start_i, start_j][0]


def init_base_to_base_score_map(score):
    base_to_map = {}
    with open(score, 'r') as f:

        reader = csv.reader(f, delimiter="\t")
        headers = reader.next()

        # Maps base to another map with scores of the bases it interacts with.
        for base in headers[1:]:
            base_to_map[base] = {}

        for row in reader:
            base = row[0]
            for i, col in enumerate(row[1:]):
                base_to_map[headers[i + 1]][base] = int(col)

    return base_to_map


# *** Main *** #
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()

    # Generator for all sequences in each fasta file.
    fasta1_generator = fastaread(command_args.seq_a)
    fasta2_generator = fastaread(command_args.seq_b)

    # We assume there is only one sequence per fasa file.
    header1, seq1 = fasta1_generator.next()
    header2, seq2 = fasta2_generator.next()

    # Initialize data structures.
    dynamic_table = np.zeros((len(seq1) + 1, len(seq2) + 1, 2))
    base_to_base_score_map = init_base_to_base_score_map(command_args.score)

    if command_args.align_type == GLOBAL:
        aligned_seq1, aligned_seq2, score = global_seq_align(dynamic_table, base_to_base_score_map, seq1, seq2, dynamic_table.shape[0], dynamic_table.shape[1])
        print_alignment(aligned_seq1, aligned_seq2, GLOBAL, score)

    elif command_args.align_type == LOCAL:
        aligned_seq1, aligned_seq2, score = local_seq_align(dynamic_table, base_to_base_score_map, seq1, seq2, dynamic_table.shape[0], dynamic_table.shape[1])
        print_alignment(aligned_seq1, aligned_seq2, LOCAL, score)

if __name__ == '__main__':
    main()
