import os
import json


CONFIG_FILE = "config.txt"



class Config:
    is_global = True
    target_list = []
    add_list = []
    ignore_list = []
    column_list = []
    max_distance = 0.9
    max_gene_amount = 100

    def __init__(self, config_raw: dict):

        self.is_global = False
        if config_raw["align_mode"] == "global":
            self.is_global = True

        self.target_list = config_raw["gene_query"]
        self.add_list = config_raw["gene_add"]
        self.ignore_list = config_raw["gene_ignore"]
        self.column_list = config_raw["columns_heatmap"]
        self.max_distance = config_raw["max_distance"]
        self.max_gene_amount = config_raw["max_gene_amount"]


def read_config():

    try:
        config_raw = json.load(open(CONFIG_FILE, 'r'))
    except json.JSONDecodeError:
        with open(CONFIG_FILE, 'w') as file:
            file.write('''
    {
        "type": "json",
        "make_list_like_this": ["Gene1", "Gene2", "Gene3"],
    
        "align_mode_explanation": "global for full sequence of protein; local for some part of it",
        "align_mode": "global",
    
        "gene_query_explanation": "Genes in the genome to find its relative genes",
        "gene_query": [],
    
        "gene_add_explanation": "Genes in the genome to add in the heatmap and tree",
        "gene_add": [],
    
        "gene_ignore_explanation": "Genes in the genome to ignore by making the heatmap and tree",
        "gene_ignore": [],
    
        "columns_heatmap_explanation": "Columns to show in the heatmap; select them from the RNA-seq file data, empty means all",
        "columns_heatmap": [],
    
        "max_distance_explanation": "Distance threshold to count the 'similar' genes, 0 means only identical / 1 means all",
        "max_distance": 0.95,
    
        "max_gene_amount_explanation": "Amount threshold of genes to show; more than 50 can make the software too slow in multiple alignment step",
        "max_gene_amount": 70,
    
        "End_of_file": true
    }
    ''')

    config_raw = json.load(open(CONFIG_FILE, 'r'))
    return Config(config_raw)






