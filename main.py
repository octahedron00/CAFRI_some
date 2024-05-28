# from project treatmap

import os
import logging
import datetime

logging.basicConfig(format='%(asctime)s / %(levelname)s : %(message)s', level=logging.INFO)

ALL_PROTEIN_FILE = 'data_genome/all_protein.fasta'
ALL_GENE_FILE = 'data_genome/all_gene.fasta'
ALL_CDS_FILE = 'data_genome/all_cds.fasta'
ALL_PROMOTER_FILE = 'data_genome/all_promoter.fasta'
ALL_PROMOTER_GENE_FILE = 'data_genome/all_promoter+gene.fasta'

DATA_GENOME = 'data_genome/'


def get_opposite_strand(seq: str):

    new_seq = ""
    for s in seq[::-1]:
        if s == 'A':
            new_seq += 'T'
        if s == 'T':
            new_seq += 'A'
        if s == 'G':
            new_seq += 'C'
        if s == 'C':
            new_seq += 'G'

    return new_seq


def read_fasta(file_name: str):

    seq_iter = SeqIO.parse(file_name, 'fasta')
    seq_map = {}

    for seq in seq_iter:
        seq_map[seq.name] = str(seq.seq).upper()

    return seq_map



def check_existing_genome_digest():

    is_good = True
    logging.info('checking genome digest files')

    if not os.path.isfile(ALL_PROTEIN_FILE):
        is_good = False
        logging.info('all_protein.fasta file not exists in data_genome')
    if not os.path.isfile(ALL_GENE_FILE):
        is_good = False
        logging.info('all_gene.fasta file not exists in data_genome')
    if not os.path.isfile(ALL_CDS_FILE):
        is_good = False
        logging.info('all_cds.fasta file not exists in data_genome')
    if not os.path.isfile(ALL_PROMOTER_FILE):
        is_good = False
        logging.info('all_promoter.fasta file not exists in data_genome')
    if not os.path.isfile(ALL_PROMOTER_GENE_FILE):
        is_good = False
        logging.info('all_promoter+gene.fasta file not exists in data_genome')

    return is_good


def write_genome_digest():
    file_genome = ''
    file_annotation = ''

    list_all = os.listdir(DATA_GENOME)
    list_files = [file for file in list_all if os.path.isfile(DATA_GENOME + file)]
    for file in list_files:
        if file.find('gff') > -1:
            file_annotation = DATA_GENOME + file
        elif (file.find('md') > -1) or (file.find('all_') > -1):
            continue
        elif file.find('fa') > -1:
            file_genome = DATA_GENOME + file

    if len(file_genome) < 1:
        logging.warning('No genome file (fasta form) was found; quitting...')
        return False
    logging.info(f"set genome file: {file_genome}")
    if len(file_annotation) < 1:
        logging.warning('No annotation file (gff form) was found; quitting...')
        return False
    logging.info(f"set annotation file: {file_annotation}")



    return True


def main():

    if not check_existing_genome_digest():
        logging.warning("genome digest files are not sufficient")
        logging.info("Run genome digest")

        if not write_genome_digest():
            return

    pass


if __name__ == '__main__':

    main()