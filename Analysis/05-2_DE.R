# differential expression with input 
library(enONE)
library(tidyverse)
library(utilsR)
library(patchwork)
library(biomaRt)

# load NAD-RNA data from results of 03_enONE.R
load('data/Enone.RData')
# get normalized counts
counts_norm <- Counts(Enone, 'sample', 'TMM_RUVse_k3')
# get raw counts
counts_sample <- Counts(Enone, 'sample', 'Raw')
# get input counts
enrich_idx <- CreateGroupMatrix(Enone$enrich)
counts_input <- counts_sample[,enrich_idx[2,]]

# perform differential expression analysis between d14.ctrl v d0.ctrl and d14.nmn v d0.nmn
# get normalized factors
enone.factor <- getFactor(Enone, 'sample', 'TMM_RUVse_k3')
deg.in <- edgeRDE(counts_sample, 
                  group = Enone$condition,
                  norm.factors = enone.factor$normFactor, # scaling factors
                  adjust.factors = enone.factor$adjustFactor, # RUV estimated factors
                  design.formula = as.formula("~0+condition"),
                  contrast.df = data.frame(Group1=c('D14_ctrl.Input','D14_nmn.Input'),
                                           Group2=c('D0_ctrl.Input','D0_nmn.Input')),
                  only.pos = FALSE)
names(deg.in$res.ls) <- c('ctrl','nmn')

# volcano plot
# differentially expressed genes is defined as absolute fold change ≥ 1.5 and 
# adjusted P < 0.05 in D14 samples compared to D0 samples.
vp1 <- utilsR::ggVolcano(deg.in$res.ls$ctrl, lfc.col='logFC', p.col='FDR', 
                         title='Water(D14/D0)', up.lfc.cutoff=log2(1.5), down.lfc.cutoff=log2(2/3))
vp2 <- utilsR::ggVolcano(deg.in$res.ls$nmn, lfc.col='logFC', p.col='FDR', 
                         title='NMN(D14/D0)', up.lfc.cutoff=log2(1.5), down.lfc.cutoff=log2(2/3)) + 
  labs(y='')

# supp. fig. 4A
vps <- vp1 + vp2 + 
  plot_annotation(title = 'Transcriptome') & 
  scale_color_manual(values=c('#124b86', '#b2b2b2', '#c77a32'))
ggsave('D14vD0_input.pdf', vps, width=7, height=4)
