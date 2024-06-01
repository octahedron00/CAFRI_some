import datetime
import math
import os

import pandas as pd
import re
import subprocess

import ete3

from Bio import Phylo, SeqIO, AlignIO
from Bio.Phylo.Consensus import bootstrap_trees, majority_consensus

from Bio.Phylo.BaseTree import *

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import PairwiseAligner, substitution_matrices

MAX_DIST = 0.9

PROTEIN_ONLY = "BCEFHIJKLMNOPQRSUVWXYZ"
BOOTSTRAP = 1

SPACER = 100
N_DATA = 1000

DIR_DATA = "data"
DIR_VARIABLE = "var"
DIR_SRC = "_src"

TASK_TITLE = str(datetime.datetime.now())[5:-10]
TASK_TITLE = "".join([c for c in TASK_TITLE if c not in "\/:*?<>| -"])


def read_fasta(file_name: str):

    seq_iter = SeqIO.parse(file_name, 'fasta')
    seq_map = {}

    for seq in seq_iter:
        seq_map[seq.name] = str(seq.seq).upper()

    return seq_map


def get_better_name(before_name: str):

    before_after_list = []

    for i in range(10):
        before_after_list.append((f"Vp{i}G", f"Vp{i:02d}G"))
        before_after_list.append((f"Vo{i}G", f"Vo{i:02d}G"))

    after_name = before_name
    for before, after in before_after_list:
        after_name = after_name.replace(before, after)

    return after_name


def get_origin_name(before_name: str):

    after_name = before_name
    after_name = after_name.replace("Vp0", "Vp")
    after_name = after_name.replace("Vo0", "Vo")
    if after_name.find(".") > 0:
        after_name = after_name[:after_name.find(".")]

    return after_name


def align_all(pos1: int, pos2: int):

    time_start = datetime.datetime.now()

    print("reading file...")
    seq_map = read_fasta('../genome/Vp_proteins.fasta')

    time_end1 = datetime.datetime.now()

    print(f"reading finished, time: {time_end1-time_start}")

    key_list = list(seq_map.keys())

    aligner = PairwiseAligner()
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -4
    aligner.left_gap_score = 0
    aligner.right_gap_score = 0
    matrix = substitution_matrices.load("BLOSUM62")

    aligner.substitution_matrix = matrix

    distance_range = [0] * (N_DATA+3)
    distance_range_for_isoforms = [0] * (N_DATA+3)
    distance_range_min = [0] * (N_DATA+3)

    time_end2 = datetime.datetime.now()
    print(f"Ready, time: {time_end2-time_end1}")

    cycle_length = int(len(seq_map.keys())/SPACER)

    total_match = int(cycle_length * (cycle_length+1) / 2)

    for i in range(cycle_length):
        n1 = (i * SPACER + pos1)
        key1 = key_list[n1]

        if len(seq_map.keys()) <= n1:
            continue

        for j in range(i, cycle_length):
            n2 = (j * SPACER + pos2)
            key2 = key_list[n2]
            if len(seq_map.keys()) <= n2:
                continue
            if n2 <= n1:
                continue

            alignments = aligner.align(seq_map[key1], seq_map[key2])

            alignment = alignments[0]

            score_max_1 = 0
            for s in seq_map[key1]:
                score_max_1 += matrix[s, s]

            score_max_2 = 0
            for s in seq_map[key2]:
                score_max_2 += matrix[s, s]

            score_min = 0

            score_alignment = alignment.score

            dist1 = (score_alignment-score_min) / (score_max_1 - score_min)
            dist1 = 1 - dist1

            dist2 = (score_alignment-score_min) / (score_max_2 - score_min)
            dist2 = 1 - dist2

            dist_int1 = int(dist1 * N_DATA)
            dist_int2 = int(dist2 * N_DATA)

            if key1[:13] == key2[:13]:
                # print(key1, key2, "as isoform")

                distance_range_for_isoforms[dist_int1] += 1
                distance_range_for_isoforms[dist_int2] += 1

            else:
                distance_range[dist_int1] += 1
                distance_range[dist_int2] += 1
                distance_range_min[min(dist_int1, dist_int2)] += 2

        delta_time = datetime.datetime.now() - time_end2

        round_now = int((i+1)*(cycle_length*2-i)/2)

        print(f"\r{round_now} / {total_match} ({round_now/total_match:.3f}), "
              f"time: {delta_time}, "
              f"estimated_time: {delta_time * ((total_match - round_now)/round_now)}         ", end="")

    time_end3 = datetime.datetime.now()
    print()
    print("File writing...")

    file_name = f"{pos1:02d}_{pos2:02d}_{TASK_TITLE}"

    with open(f"dist_range_{file_name}.txt", 'w') as file:
        file.write(f"{file_name}\n")
        for i in range(N_DATA+1):
            file.write(f"{distance_range[i]}\n")

    with open(f"dist_range_isoform_{file_name}.txt", 'w') as file:
        file.write(f"{file_name}\n")
        for i in range(N_DATA+1):
            file.write(f"{distance_range_for_isoforms[i]}\n")

    with open(f"dist_range_min_{file_name}.txt", 'w') as file:
        file.write(f"{file_name}\n")
        for i in range(N_DATA+1):
            file.write(f"{distance_range_min[i]}\n")

    print(f"File written, time: {datetime.datetime.now() - time_end3}, total: {datetime.datetime.now() - time_start}")

    return file_name, distance_range, distance_range_min, distance_range_for_isoforms


if __name__ == '__main__':
    distance_range_list = []
    distance_range_min_list = []
    distance_range_for_isoforms_list = []
    file_name_list = []
    for n in range(20):
        fn, dr, drm, dri = align_all(20+n, 90)

        file_name_list.append(fn)
        distance_range_list.append(dr)
        distance_range_min_list.append(drm)
        distance_range_for_isoforms_list.append(dri)

    with open(f"total_dist_range_{TASK_TITLE}.txt", 'w') as file:
        for i in range(len(file_name_list)):
            file.write(f"{file_name_list[i]}\t")
        file.write("\n")
        for n in range(N_DATA + 1):
            for i in range(len(file_name_list)):
                file.write(f"{distance_range_list[i][n]}\t")
            file.write("\n")

    with open(f"total_dist_range_isoform_{TASK_TITLE}.txt", 'w') as file:
        for i in range(len(file_name_list)):
            file.write(f"{file_name_list[i]}\t")
        file.write("\n")
        for n in range(N_DATA + 1):
            for i in range(len(file_name_list)):
                file.write(f"{distance_range_for_isoforms_list[i][n]}\t")
            file.write("\n")

    with open(f"total_dist_range_min_{TASK_TITLE}.txt", 'w') as file:
        for i in range(len(file_name_list)):
            file.write(f"{file_name_list[i]}\t")
        file.write("\n")
        for n in range(N_DATA + 1):
            for i in range(len(file_name_list)):
                file.write(f"{distance_range_min_list[i][n]}\t")
            file.write("\n")



