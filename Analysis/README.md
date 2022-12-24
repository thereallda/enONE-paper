## 01_RNA-seq_pipeline

In "01-RNA-seq_pipeline" folder:

- `RNAseq_pipeline.sh`: In-house script for RNA-sequencing data quality control, mapping and quantification. 
- `README.md`: Detailed description of the `RNAseq_pipeline.sh`



## 02_SequencingStats

`02_SequencingStats.R`: for assessment of sequencing statistics, including alignment summary, sequencing saturation, and synthetics RNA enrichment based on raw counts. 



## 03_enONE

`03_enONE.R`:  using enONE for data normalization and NAD-RNA identification. 



## 04_NAD-RNA_characterization

`04_NAD-RNA_characterization.R`: for NAD-RNA characterization, including gene type, chromosome distribution, and gene length. 



## 05-1_NAD-RNA_dynamics

`05-1_NAD-RNA_dynamics.R`: for differential NAD-RNA modification analysis. 



## 05-2_DE

`05-2_DE.R`: for differential gene expression analysis. 



## 06_NAD-RNA_biomarkers

`06_NAD-RNA_biomarkers.R`: for identification of NAD-RNA biomarkers. We construct a logistic regression model with an elastic net penalty based on the scaled NAD modification levels. 