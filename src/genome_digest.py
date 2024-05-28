
import os
import logging
import csv
import datetime

from Bio import SeqIO, Seq


logging.basicConfig(format='%(asctime)s / %(levelname)s : %(message)s', level=logging.INFO)


ALL_PROTEIN_FILE = 'data_genome/all_protein.fasta'
ALL_GENE_FILE = 'data_genome/all_gene.fasta'
ALL_CDS_FILE = 'data_genome/all_cds.fasta'
ALL_PROMOTER_FILE = 'data_genome/all_promoter.fasta'

DATA_GENOME = 'data_genome/'

PROMOTER_LENGTH = 2000


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


def write_fasta(file_name: str, seq_map: dict):

    with open(file_name, 'w') as file:
        for key, item in seq_map.items():
            file.write(f">{key}\n{item}\n")


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

    return is_good


def find_or_max(position):
    if position < 0:
        return 100000
    return position


def get_sequence(seq: str, strand: str):
    if strand == '-':
        return get_opposite_strand(seq)
    return seq


def try_genome_digest():
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

    '''
    now, start making files
    '''

    do_genome_digest(file_genome, file_annotation)
    return True


def do_genome_digest(file_genome, file_annotation):

    annotation_reader = csv.reader(open(file_annotation, 'r'), delimiter='\t')
    genome = read_fasta(file_genome)

    dict_gene = dict()
    dict_cds = dict()
    dict_protein = dict()
    dict_promoter = dict()

    for row in annotation_reader:
        if row[0][0] == '#' or len(row) < 9:
            print(row)
            continue
        seq_id = row[0]
        seq_type = row[2]
        start = int(row[3]) - 1
        end = int(row[4])
        # start:end < will match for all strands
        strand = row[6]
        desc = row[8]

        name = desc[desc.find("ID=")+3 : min(find_or_max(desc.find(';')),
                                             find_or_max(desc.find('.exon')),
                                             find_or_max(desc.find('.CDS')))]

        # if seq_type == 'start_codon':
        #     print(get_sequence(genome[seq_id][start:end], strand=strand))

        sequence = get_sequence(genome[seq_id][start:end], strand=strand)

        if seq_type == 'gene':
            dict_gene[name] = sequence
            if strand == '+':
                promoter = get_sequence(genome[seq_id][max(0, start-PROMOTER_LENGTH):start], strand=strand)
            else:
                promoter = get_sequence(genome[seq_id][end:min(end+PROMOTER_LENGTH, len(genome[seq_id]))], strand=strand)
            dict_promoter[name] = promoter

        if seq_type == 'CDS':
            if name not in dict_cds.keys():
                dict_cds[name] = sequence
            elif strand == '+':
                dict_cds[name] = dict_cds[name] + sequence
            elif strand == '-':
                dict_cds[name] = sequence + dict_cds[name]

    for name, sequence in dict_cds.items():
        dict_protein[name] = Seq.translate(sequence)

    write_fasta(ALL_PROTEIN_FILE, dict_protein)
    write_fasta(ALL_GENE_FILE, dict_gene)
    write_fasta(ALL_CDS_FILE, dict_cds)
    write_fasta(ALL_PROMOTER_FILE, dict_promoter)

