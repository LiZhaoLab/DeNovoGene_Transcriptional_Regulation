# Transcriptional regulation of de novo genes in D. melanogaster

This repo contains files and scripts to:  
1. Filter and extract expression matrix and compute z-score from fly cell atlas dataset.
2. Construct the transcriptional regulation network for de novo genes in D. melanogaster.
  
## Scripts  
`combined_regulon_motif.py`
- combine individual Regulon analysis into a final one.
  1. The master transcription factor must appear at least 80% of the predictions to ensure a high level of consensus.  
  2. The target genes of the master TFs must appear at least five times to make sure that they do not appear randomly.  

`extract_fca_expression.py`  
- extract sub expression matrix that:  
  1. includes only the counts of transcription factors, de novo gene candidates, and testis-biased genes in each cell, and   
  2. filtered out cells where no de novo genes were expressed and transcription factors were lowly expressed.  
- analyze and plot the cell numbers and ratios of original and final expression matrix  
- compute zscore matrix by different features: `tissue`,`annotation_broad`, and `annotation`. The three features corresponds to different tissue types, major cell types, and more specific cell types.

`network_analysis.py`  
- visualize the transcriptional regulatory network.  
- analyze the complexity of transcriptional regulations. 

`regression_and_AIC_analysis.R`  
- regression analysis between vis expression and other target genes across different RNA-seq samples.  
- AIC model for different regression analysis. the AIC analysis identified T~V as the best fit model.  
- visulization of the regression analysis. the figures including Figure 4 and Figure S8 in our original manuscript
  
## Required files  
   ```
   all_fca_genes.validation.csv  
   fca_testis_biased_genes.csv
   denovo_candidates.name.csv
   denovo_candidates.name.txt
   denovo_regulons_combined.pickle
   raw_predictions/
   ``` 
  
Note: The folder `raw_predictions` contains the individual pyscenic Regulon analysis.  
Additional datasets that are required can be downloaded from other resources:  
1. `r_fca_biohub_all_wo_blood_10x.loom` can be downloaded from Fly Cell Atlas (https://flycellatlas.org) and
2. `allTFs_dmel.txt` can be downloaded from cisTarget database (https://resources.aertslab.org/cistarget/databases/)  
3. We have included some essential files from cisTarget database in this repo, including `allTFs_dmel.txt`, `dm6_v10_clust.genes_vs_motifs.rankings.feather`, and `dm6_v10_clust.genes_vs_motifs.scores.feather`

## Required packages  
python  
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
R  
```
Rstudio=2024.09.0+375
R=4.4.1
ggplot2
ggpubr
car
dplyr
plyr
data.table
stringr 
reshape2
qvalue
scales
ggiraph
ggiraphExtra
plotrix   
AICcmodavg
cowplot
```
