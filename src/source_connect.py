import logging
import datetime
import math
import os
import subprocess

import numpy
import pandas as pd

from src.protein_compare import Protein
from Bio import Phylo, SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from src.blosum62.protein_matrices import protein_matrices, DistanceCalculator


DATAFRAME_FILE = "temp/_df.csv"
METADATA_FILE = "temp/_md.csv"

BLOSUM62 = "src/blosum62/blosum62.txt"

DIR_SRC = "src"
DIR_VARIABLE = "temp"

IS_LOG = True


def get_dataframe_all(file_list: list[str], column_list: list[str]):
    if len(file_list) < 1:
        return pd.DataFrame([[1, 2], [1, 2]], index=["a", "b"], columns=["null_1", "null_2"])

    df = pd.read_csv(file_list[0], index_col=0)
    logging.info("read_1 complete")

    for i, file_address in enumerate(file_list[1:]):
        df2 = pd.read_csv(file_address, index_col=0)
        df = pd.concat([df, df2], axis=1)
        logging.info(f"read_{i + 2} complete + concatenated")

    if len(column_list) > 0:
        df = df.loc[:, column_list]

    if len(df.columns) < 2:
        df2 = pd.DataFrame([1, 2], index=["a", "b"], columns=["null"])
        df = pd.concat([df, df2], axis=1)

    # df = df.fillna(0)

    # print(df)
    return df


def make_dataframe_file(valid_protein_list: list[Protein], file_address_list: list[str], column_list: list[str]):
    df = get_dataframe_all(file_address_list, column_list)

    n_column = len(df.columns)

    with open(DATAFRAME_FILE, 'w') as file:

        file.write("Gene")
        for column_name in df.columns.values:
            file.write(f",{column_name}")
        file.write("\n")

        for protein in valid_protein_list:
            if protein.gene_name not in df.index:
                file.write(protein.get_best_name())
                for i in range(n_column):
                    file.write(f",{numpy.NaN}")
                file.write("\n")
                continue

            row = df.loc[protein.gene_name]
            # print(row)

            file.write(protein.get_best_name())
            for column_name in df.columns.values:

                if IS_LOG:
                    file.write(f",{math.log10(row[column_name] + 1)}")
                else:
                    file.write(f",{row[column_name]}")
            file.write("\n")


def make_metadata_file(valid_protein_list: list[Protein], target_list: list[str]):
    with open(METADATA_FILE, 'w') as file:

        file.write(f"name,group\n")

        for protein in valid_protein_list:
            file.write(f"{protein.get_best_name()},{protein.metadata}\n")


def align_by_clustalw():

    command = f"{os.path.join(DIR_SRC, 'clustalw2', 'clustalw2.exe')} " \
              f"-infile=\"{os.path.join(DIR_VARIABLE, '_align.fasta')}\" " \
              f"-matrix=BLOSUM"

    logging.info(f"command: {command}")
    result = subprocess.run(command, shell=True)
    logging.info(f"error: {result.stderr}")


def make_tree():
    multiple_align = AlignIO.read(os.path.join(DIR_VARIABLE, '_align.aln'), "clustal")
    AlignIO.write(multiple_align, os.path.join(DIR_VARIABLE, '_align.phy'), 'phylip-relaxed')

    calculator = DistanceCalculator(model_address=BLOSUM62)
    constructor = DistanceTreeConstructor()

    distance_matrix = calculator.get_distance(multiple_align)

    tree = constructor.nj(distance_matrix)

    print(tree)

    Phylo.write(tree, os.path.join(DIR_VARIABLE, "tree.xml"), 'phyloxml')
    Phylo.write(tree, os.path.join(DIR_VARIABLE, "_tree.dnd"), 'newick')

    # command = f"{os.path.join(DIR_SRC, 'fasttree', 'FastTree.exe')} -mlacc 2 -slownni " \
    #           f"{os.path.join(DIR_VARIABLE, '_align.phy')} > {os.path.join(DIR_VARIABLE, '_tree.dnd')}"

    # logging.info(f"command: {command}")
    # result = subprocess.run(command, shell=True)
    # logging.info(f"error: {result.stderr}")


def make_treatmap_image_by_r():
    command = f"{os.path.join(DIR_SRC, 'compact_R', 'bin', 'Rscript.exe')} " \
              f"{os.path.join(DIR_SRC, 'R_script', 'draw_tree.R')}"

    logging.info(f"command: {command}")
    result = subprocess.run(command, shell=True)
    logging.info(f"error: {result.stderr}")


if __name__ == '__main__':
    DIR_SRC = "..\\src"
    DIR_VARIABLE = "..\\temp"
    align_by_clustalw()
    make_tree()
    make_treatmap_image_by_r()
