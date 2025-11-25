# Heimdall_CCA_paper
pipelines and scripts used in the workflow 

# motif_PAP_vs_CCA.py
input fasta file can contain all CCA sequences (ie a concatenated file for all), and returns a csv file with headers:
- Protein: this is the accession of the given sequence 
- Enzyme Type: if it classified as "CCA, PAP, or Unknown"
- Missing Motifs: lists all/any missing motifs 
- Insufficient Upstream Space: defines if there was insufficient upstream space to identify the characteristic PAP upstream motif
- Final Classification: Given a set of conditions, was it possible to confidently assign the proteins/enzyme to a specific type
