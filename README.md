# Transcriptional regulation of de novo genes in D. melanogaster

This repo contains files and scripts to:  
1. Filter and extract expression matrix from fly cell atlas dataset.
2. Construct the transcriptional regulation network for de novo genes in D. melanogaster.
  
Scripts  
`extract_fca_expression.py`  
- extract sub expression matrix that (a) includes only the counts of transcription factors, de novo gene candidates, and testis-biased genes in each cell, and (b) filtered out cells where no de novo genes were expressed and transcription factors were lowly expressed.  
- analyze and plot the cell numbers and ratios of original and final expression matrix  

`network_analysis.py`  
- visualize the transcriptional regulatory network.  
- analyze the complexity of transcriptional regulations.  
  
Required files  
    `all_fca_genes.validation.csv`  
    `fca_testis_biased_genes.csv`  
    `denovo_candidates.name.csv`  
    `denovo_candidates.name.txt`  
    `denovo_regulons_combined.pickle`  
  
Note: Additional dataset can be downloaded from Fly Cell Atlas (https://flycellatlas.org) and cisTarget database (https://resources.aertslab.org/cistarget/databases/)

Required packages  
```
python=3.7.12
numpy=1.19.5
seaborn=0.12.2
networkx=2.5
matplotlib=3.5.3
scanpy=1.9.1
scipy=1.7.1
pandas=1.3.5
```
