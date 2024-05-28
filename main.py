# from project treatmap

import os
import logging
import csv
import datetime

from Bio import SeqIO, Seq

from src.genome_digest import check_existing_genome_digest, try_genome_digest, read_fasta
from src.protein_compare import get_query_protein_list, get_total_protein_list, get_similar_protein_list

logging.basicConfig(format='%(asctime)s / %(levelname)s : %(message)s', level=logging.INFO)

ALL_PROTEIN_FILE = 'data_genome/all_protein.fasta'
ALL_GENE_FILE = 'data_genome/all_gene.fasta'
ALL_CDS_FILE = 'data_genome/all_cds.fasta'
ALL_PROMOTER_FILE = 'data_genome/all_promoter.fasta'

QUERY_PROTEIN_FILE = 'query.txt'

MIN_ROW = 1
MAX_ROW = 500


def main():
    if not check_existing_genome_digest():
        logging.warning("genome digest files are not sufficient")
        logging.info("Run genome digest")

        if not try_genome_digest():
            return

    target_list = []
    add_list = []
    ignore_list = []
    column_list = []
    rna_seq_file_list = []

    logging.info(f"target {target_list}")
    logging.info(f"add {add_list}")
    logging.info(f"ignore {ignore_list}")

    all_protein_dict = read_fasta(ALL_PROTEIN_FILE)
    all_gene_dict = read_fasta(ALL_GENE_FILE)
    all_cds_dict = read_fasta(ALL_CDS_FILE)
    all_promoter_dict = read_fasta(ALL_PROMOTER_FILE)

    query_protein_dict = read_fasta(QUERY_PROTEIN_FILE)

    target_protein_list = get_query_protein_list(target_list, all_protein_dict, query_protein_dict)
    total_protein_list = get_total_protein_list(all_protein_dict, all_gene_dict, all_cds_dict, all_promoter_dict,
                                                query_protein_dict)

    sorted_protein_list = get_similar_protein_list(target_protein_list, total_protein_list, add_list, ignore_list)
    logging.info(f"{len(sorted_protein_list)} found")


if __name__ == '__main__':

    main()

