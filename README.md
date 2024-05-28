# CAFRI-some

CRISPR Applicable Functional Redundancy Inspector

to accelerate functional genomics in some/any genome

Similar scheme in [CAFRI-Rice](https://cafri-rice.khu.ac.kr/)

## Work scheme

1. User inputs total genome(fasta) & annotation file(gff), and RNA-seq result file(?, optional)
2. software makes all_protein, gene, cds, promoter, gene-promoter fasta file
3. User inputs query of...

- Any protein sequence to find the similar proteins in genome
- Any gene name list to find the similar proteins in genome
- Any gene name to add / ignore from the result of similar proteins

4. Software finds all the similar proteins throughout the whole genome 
5. Software draws the phylogenetic tree(maximum_likelihood) from the protein sequences
6. Software shows the tree (with a heatmap of the RNA-seq result, optional)

## Dependencies

- SeqIO: for reading fasta files
