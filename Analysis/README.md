## 01_RNA-seq_pipeline

In "01-RNA-seq_pipeline" folder:

- `RNAseq_pipeline.sh`: In-house script for RNA-sequencing data quality control, mapping and quantification. 
- `README.md`: Detailed description of the `RNAseq_pipeline.sh`



## 02_SequencingStats

`02_SequencingStats.R`: for assessment of sequencing statistics, including alignment summary, and sequencing saturation.



## 03_enONE

`03_enONE.R`:  using enONE for data normalization and NAD-RNA identification. 



## 04_NAD-RNA_characterization

`04-1_NAD-RNA_characterization.R`: for NAD-RNA characterization, including gene type, chromosome distribution, gene length, and intron numbers. 

`04-2_NAD-RNA_MFE.R`: for analysis of minimum free energy (MFE) based on 5'-UTR sequence in R

`04-2_RNAFold.sh`: used `RNAfold` to calculate the minimum free energy (MFE) based on 5'-UTR sequence in shell.

`04-2_RNAFold.sh` can be used in the following command 

> Prerequisite: ViennaRNA package (https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html)
>
> Before using `04-2_RNAFold.sh`, one should get sequences of interest by genome coordinations (bed format), which can be done with the following commands
>
> ```bash
> bedtools getfasta -name+ -fi GRCh38.primary_assembly.genome.fa -bed NADRNA_5UTR.bed -fo NADRNA_5UTR.fa
> ```
>
> `GRCh38.primary_assembly.genome.fa` should be indexed and the index should be placed in the same directory as the genome fasta file. 

```bash
bash 04-2_RNAFold.sh -i NADRNA_5UTR.fa -o NADRNA_MFE.csv
```

`-i`: input sequence extracted based on the 5'UTR coordination, in fasta format

`-o`: MFE of each sequence. `04-2_NAD-RNA_MFE.R` use this file for visualization



## 05_NAD-RNA_dynamics

`05_NAD-RNA_dynamics.R`: for assessing NAD-RNA dynamics with age. 



## 06_AgingClock

`06_AgingClock.R`: for building age prediction models (elastic net, lasso, ridge regression models) and validation of age prediction model that combined signatures from transcriptome and epitranscriptome. 



`helper.R`: custom functions for analysis.