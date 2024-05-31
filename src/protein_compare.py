import datetime
import logging
import math
import os

from Bio.Align import PairwiseAligner, substitution_matrices

COMPARED_PROTEIN_FILE = "result/compared_protein.fasta"
COMPARED_GENE_FILE = "result/compared_gene.fasta"
COMPARED_CDS_FILE = "result/compared_cds.fasta"
COMPARED_PROMOTER_FILE = "result/compared_promoter.fasta"
COMPARED_PROMOTER_GENE_FILE = "result/compared_promoter+gene.fasta"
COMPARED_RESULT_FILE = "result/compared_result.fasta"
FILE_FOR_MULTIPLE_ALIGN = "temp/_align.fasta"


def write_fasta(file_name: str, seq_map: dict):

    with open(file_name, 'w') as file:
        for key, item in seq_map.items():
            file.write(f">{key}\n{item}\n")


class Protein:
    name = ""
    gene_name = ""
    seq = ""

    gene_seq = ""
    cds_seq = ""
    promoter_seq = ""

    distance_main = 10.0
    align_result = ""
    align_score = ""
    align_score_base = ""
    align_with = ""
    metadata = ""
    is_isoform = False

    def __init__(self, key: str, seq: str, is_external=False):

        self.name = key
        self.seq = seq.upper().replace('B', '').replace('J', '').replace('O', '').replace('U', '').replace('X', '').replace('Y', '').replace('Z', '')
        self.gene_seq = ""
        self.cds_seq = ""
        self.promoter_seq = ""

        if is_external:
            self.gene_name = key
        else:
            self.gene_name = key[:key.find('.')]
        self.distance = 10.0
        self.metadata = "none"
        self.is_isoform = False

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return f">{self.get_best_name()} ({self.distance:6f})\n" \
               f"{self.seq}\n"

    def __lt__(self, other):
        return self.distance < other.distance

    def set_distance(self, distance: float):
        self.distance = distance

    def set_align_result(self, align_result: str, align_score: int, align_score_base: int, align_with: str):
        self.align_result = align_result
        self.align_score = align_score
        self.align_score_base = align_score_base
        self.align_with = align_with

    def set_isoform(self):
        self.is_isoform = True

    def set_metadata(self, metadata: str):
        self.metadata = metadata

    def set_gene_cds_promoter(self, gene_seq: str, cds_seq: str, promoter_seq: str):
        self.gene_seq = gene_seq
        self.cds_seq = cds_seq
        self.promoter_seq = promoter_seq

    def get_best_name(self):
        if self.is_isoform:
            return self.name
        else:
            return self.gene_name


def get_query_protein_list(target_name_list: list[str], all_protein_dict: dict, query_protein_dict: dict):
    query_protein_list = []

    for key in all_protein_dict.keys():
        for target_name in target_name_list:
            if target_name.find('.') > -1 and key == target_name:
                protein = Protein(key=key, seq=all_protein_dict[key])
                query_protein_list.append(protein)
            elif key[:key.find('.')] == target_name:
                protein = Protein(key=key, seq=all_protein_dict[key])
                query_protein_list.append(protein)

    for key in query_protein_dict.keys():
        print(key)
        ext_protein = Protein(key=key, seq=query_protein_dict[key])
        query_protein_list.append(ext_protein)

    return get_query_protein_list_with_distance(query_protein_list)


def get_query_protein_list_with_distance(query_protein_list: list[Protein]):

    query_protein_list_with_distance = []

    aligner = PairwiseAligner()
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -4
    aligner.left_gap_score = 0
    aligner.right_gap_score = 0
    matrix = substitution_matrices.load("BLOSUM62")

    for protein_1 in query_protein_list:
        protein_1.set_distance(0)
        score_max = get_aligned_score(aligner, matrix, protein_1, protein_1, is_self=True)
        score_min = 0

        for protein_2 in query_protein_list:
            if protein_1.seq == protein_2.seq:
                continue

            score_alignment, _ = get_aligned_score(aligner, matrix, protein_1, protein_2)
            distance = 1 - (score_alignment - score_min) / (score_max - score_min)

            if distance > protein_1.distance:
                protein_1.set_distance(distance)
                print(distance)

        query_protein_list_with_distance.append(protein_1)
    return query_protein_list_with_distance


def get_total_protein_list(all_protein_dict: dict, all_gene_dict: dict, all_cds_dict: dict, all_promoter_dict: dict,
                           query_protein_dict: dict):
    total_protein_list = []

    for key in sorted(list(all_protein_dict.keys())):
        protein = Protein(key=key, seq=all_protein_dict[key])

        if int(key[key.find(".") + 1:]) > 1:
            protein.set_isoform()

            previous_protein = total_protein_list[-1]
            previous_protein.set_isoform()
            total_protein_list = total_protein_list[:-1] + [previous_protein]

        # print(protein.name, protein.gene_name)
        protein.set_gene_cds_promoter(all_gene_dict[protein.gene_name], all_cds_dict[protein.name], all_promoter_dict[protein.gene_name])
        total_protein_list.append(protein)

    for key in query_protein_dict.keys():
        ext_protein = Protein(key=key, seq=query_protein_dict[key], is_external=True)
        total_protein_list.append(ext_protein)

    return total_protein_list


def get_aligned_score(aligner: PairwiseAligner, matrix, query_protein: Protein, target_protein: Protein, is_self=False):

    aligner.substitution_matrix = matrix

    if is_self:
        score = 0
        for s in query_protein.seq:
            score += matrix[s, s]
        return score

    # print(query_protein.seq, target_protein.seq)
    alignments = aligner.align(query_protein.seq, target_protein.seq)

    alignment = alignments[0]
    match_line = "" # nottodo
    alignment_str = f"{alignment[0]}\n{match_line}\n{alignment[1]}"

    return alignment.score, alignment_str


def get_similar_protein_list(query_protein_list: list[Protein], total_protein_list: list[Protein],
                             add_list: list[str], ignore_list: list[str], max_distance=0.9, max_gene_amount=100):

    MAX_DIST = max_distance
    MAX_GENE_AMOUNT = max_gene_amount

    time_start = datetime.datetime.now()

    aligner = PairwiseAligner()
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -4
    aligner.left_gap_score = 0
    aligner.right_gap_score = 0
    matrix = substitution_matrices.load("BLOSUM62")

    query_protein_names_list = [query_protein.name for query_protein in query_protein_list]

    for i, total_protein in enumerate(total_protein_list):

        if total_protein.name in ignore_list or total_protein.gene_name in ignore_list:
            total_protein.set_distance(10)
            total_protein.set_align_result(align_result="ignored_manually\n", align_score=0,
                                           align_score_base=0, align_with="none")
            total_protein_list[i] = total_protein
            continue
        # if len(query_protein_list) < 1:
        #     continue

        score_max_gene = get_aligned_score(aligner, matrix, total_protein, total_protein, is_self=True)

        for j, query_protein in enumerate(query_protein_list):

            if j < 1 and total_protein.name in query_protein_names_list:
                total_protein.set_metadata("origin")
                total_protein.set_distance(0)
                total_protein.set_align_result(total_protein.seq + "\n" + "" + "\n" + total_protein.seq + "\n",
                                               score_max_gene, score_max_gene, total_protein.name)
                break

            score_min = 0 * min(len(total_protein), len(query_protein))

            score_max_query = get_aligned_score(aligner, matrix, query_protein, query_protein, is_self=True)
            score_alignment, alignment_str = get_aligned_score(aligner, matrix, query_protein, total_protein)

            distance = 1 - (score_alignment - score_min) / (score_max_query - score_min)

            if total_protein.distance > distance:
                if distance < MAX_DIST:
                    total_protein.metadata = "align"
                total_protein.set_distance(distance=distance)
                total_protein.set_align_result(align_result=alignment_str, align_score=score_alignment,
                                               align_score_base=score_max_query, align_with=query_protein.name)

            if distance > (query_protein.distance + MAX_DIST):
                break

        if total_protein.distance < 0.2:
            print()
            logging.info(f"Found identical-ish: {total_protein.name} ({total_protein.distance}; {total_protein.seq})")

        if total_protein.name in add_list or total_protein.gene_name in add_list:
            # print("add", total_protein.name, total_protein.gene_name)
            if total_protein.metadata.find('align') > -1:
                total_protein.metadata = "align+add"
            else:
                total_protein.metadata = "add"
            total_protein.set_distance(-1)
            total_protein.set_align_result(align_result="added_manually\n", align_score=0,
                                           align_score_base=0, align_with="none")

        total_protein_list[i] = total_protein
        time_now = datetime.datetime.now()
        print(f"\r{i + 1}/{len(total_protein_list)} : {total_protein.distance:.03f} "
              f"({time_now - time_start} passed / {((time_now - time_start) / (i + 1)) * (len(total_protein_list) - i - 1)} left)",
              end="")

    sorted_protein_list = sorted(total_protein_list)

    sorted_protein_list = [protein for protein in sorted_protein_list if protein.distance < MAX_DIST]
    sorted_protein_list = sorted_protein_list[:min(len(sorted_protein_list), MAX_GENE_AMOUNT)]

    compared_protein_dict = {f"{protein.name} {protein.metadata} ({protein.distance:.6f})": protein.seq
                             for protein in sorted_protein_list}
    compared_gene_dict = {f"{protein.name} {protein.metadata} ({protein.distance:.6f})": protein.gene_seq
                          for protein in sorted_protein_list}
    compared_cds_dict = {f"{protein.name} {protein.metadata} ({protein.distance:.6f})": protein.cds_seq
                         for protein in sorted_protein_list}
    compared_promoter_dict = {f"{protein.name} {protein.metadata} ({protein.distance:.6f})": protein.promoter_seq
                              for protein in sorted_protein_list}
    compared_promoter_gene_dict = {f"{protein.name} {protein.metadata} ({protein.distance:.6f})": protein.promoter_seq + protein.gene_seq
                                   for protein in sorted_protein_list}
    compared_result_dict = {f"{protein.name} {protein.metadata} ({protein.distance:.6f}) with {protein.align_with}":
                            protein.align_result for protein in sorted_protein_list}

    compared_protein_for_multiple_align_dict = {protein.get_best_name(): protein.seq for protein in sorted_protein_list}

    write_fasta(FILE_FOR_MULTIPLE_ALIGN, compared_protein_for_multiple_align_dict)
    write_fasta(COMPARED_PROTEIN_FILE, compared_protein_dict)
    write_fasta(COMPARED_GENE_FILE, compared_gene_dict)
    write_fasta(COMPARED_CDS_FILE, compared_cds_dict)
    write_fasta(COMPARED_PROMOTER_FILE, compared_promoter_dict)
    write_fasta(COMPARED_PROMOTER_GENE_FILE, compared_promoter_gene_dict)
    write_fasta(COMPARED_RESULT_FILE, compared_result_dict)

    return sorted_protein_list
