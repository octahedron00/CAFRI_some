# from project treatmap

import os
import logging
import csv
import datetime
import shutil

from Bio import SeqIO, Seq

from _src.genome_digest import check_existing_genome_digest, try_genome_digest, read_fasta
from _src.protein_compare import get_query_protein_list, get_total_protein_list, get_similar_protein_list
from _src.source_connect import make_dataframe_file, make_metadata_file, align_by_clustalw, make_tree, make_treatmap_image_by_r
from _src.config import Config, read_config

logging.basicConfig(format='%(asctime)s / %(levelname)s : %(message)s', level=logging.INFO)

ALL_PROTEIN_FILE = 'data_genome/all_protein.fasta'
ALL_GENE_FILE = 'data_genome/all_gene.fasta'
ALL_CDS_FILE = 'data_genome/all_cds.fasta'
ALL_PROMOTER_FILE = 'data_genome/all_promoter.fasta'

QUERY_PROTEIN_FILE = 'query.txt'
CONFIG_FILE = "config.txt"
ALL_RNA_SEQ_FILE_FOLDER = 'data_rna_seq'

VERSION = "1.12 b.2024.06.03"

TASK_TITLE = "task " + str(datetime.datetime.now())[5:-10]


def main():
    if not check_existing_genome_digest():
        logging.warning("genome digest files are not sufficient")
        logging.info("Run genome digest")

        if not try_genome_digest():
            return

    result_folder_name = "result_" + "".join([c for c in TASK_TITLE if c not in "\/:*?<>| -"])
    config = read_config()
    rna_seq_file_list = [os.path.join(ALL_RNA_SEQ_FILE_FOLDER, file) for file in os.listdir(ALL_RNA_SEQ_FILE_FOLDER)
                         if os.path.isfile(os.path.join(ALL_RNA_SEQ_FILE_FOLDER, file))]

    os.mkdir(result_folder_name)
    shutil.copyfile(CONFIG_FILE, os.path.join(result_folder_name, CONFIG_FILE))
    shutil.copyfile(QUERY_PROTEIN_FILE, os.path.join(result_folder_name, QUERY_PROTEIN_FILE))

    logging.info(f"target {config.target_list}")
    logging.info(f"add {config.add_list}")
    logging.info(f"ignore {config.ignore_list}")
    logging.info(f"rna_seq files: {rna_seq_file_list}")
    logging.info(f"columns selected: {config.column_list}")
    logging.info(f"global align: {config.is_global}")

    all_protein_dict = read_fasta(ALL_PROTEIN_FILE)
    all_gene_dict = read_fasta(ALL_GENE_FILE)
    all_cds_dict = read_fasta(ALL_CDS_FILE)
    all_promoter_dict = read_fasta(ALL_PROMOTER_FILE)

    query_protein_dict = read_fasta(QUERY_PROTEIN_FILE)

    target_protein_list = get_query_protein_list(config.target_list, all_protein_dict, query_protein_dict,
                                                 is_global=config.is_global)
    total_protein_list = get_total_protein_list(all_protein_dict, all_gene_dict, all_cds_dict, all_promoter_dict,
                                                query_protein_dict)

    sorted_protein_list = get_similar_protein_list(target_protein_list, total_protein_list, config.add_list,
                                                   config.ignore_list, config.max_distance, config.max_gene_amount,
                                                   is_global=config.is_global, folder=result_folder_name)
    logging.info(f"finished")
    logging.info(f"{len(sorted_protein_list)} found")

    if len(sorted_protein_list) < 2:
        sorted_protein_list = sorted_protein_list + sorted_protein_list

    make_dataframe_file(sorted_protein_list, rna_seq_file_list, config.column_list)
    make_metadata_file(sorted_protein_list, config.target_list)

    align_by_clustalw()

    make_tree()

    make_treatmap_image_by_r(result_folder_name)


if __name__ == '__main__':
    # make_tree()
    # make_treatmap_image_by_r()

    print(f"CAFRI-some ver. {VERSION}")
    main()

