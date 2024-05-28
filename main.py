# from project treatmap

import os
import logging
import csv
import datetime

from Bio import SeqIO, Seq

from src.genome_digest import check_existing_genome_digest, write_genome_digest

logging.basicConfig(format='%(asctime)s / %(levelname)s : %(message)s', level=logging.INFO)

ALL_PROTEIN_FILE = 'data_genome/all_protein.fasta'
ALL_GENE_FILE = 'data_genome/all_gene.fasta'
ALL_CDS_FILE = 'data_genome/all_cds.fasta'
ALL_PROMOTER_FILE = 'data_genome/all_promoter.fasta'
ALL_PROMOTER_GENE_FILE = 'data_genome/all_promoter+gene.fasta'


def main():
    if not check_existing_genome_digest():
        logging.warning("genome digest files are not sufficient")
        logging.info("Run genome digest")

        if not write_genome_digest():
            return

    pass


if __name__ == '__main__':

    main()