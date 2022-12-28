# enONE for data normalization and NAD-RNA identification
library(enONE)
library(tidyverse)
library(patchwork)

# In Data folder 
metadata <- read.csv('Data/metadata.csv', comment.char = '#')
counts_df <- read.csv('Data/Counts.csv', row.names = 1)

# prefix of Drosophila spike-in gene id
spikeInPrefix <- '^FB'
# create group matrix for RUV
sc_idx <- CreateGroupMatrix(metadata$condition)
enrich_group <- str_extract(metadata$condition, "(Input)|(Enrich)") # get id for enrichment group
enrich_idx <- CreateGroupMatrix(enrich_group)

# filtering ----
## filter low-expressed genes
## kept genes with at least 0.44 CPM in at least 4 samples
keep <- edgeR::filterByExpr(counts_df, group = metadata$condition, min.count = 20)
counts_keep <- counts_df[keep,]

## filter rRNA of dog
## annotation retrieved from https://ftp.ensembl.org/pub/release-106/gff3/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz
gtf.cf <- rtracklayer::import('Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf')
gene.cf <- as.data.frame(subset(gtf.cf, type == 'gene'))
rrna.id <- subset(gene.cf, gene_biotype == 'rRNA')$gene_id
counts_keep <- counts_keep[!rownames(counts_keep) %in% rrna.id, ]

## ronser's test for outlier assessment
vt <- DESeq2::vst(as.matrix(counts_keep), nsub = 20000)
pc <- prcomp(t(vt))
## all samples pass the test
EnvStats::rosnerTest(pc$x[,1])

# create Enone ----
Enone <- createEnone(counts_keep, 
                     bio.group = metadata$condition,
                     enrich.group = enrich_group,
                     batch.group = NULL,
                     spike.in.prefix = spikeInPrefix,
                     synthetic.id = c('Syn1','Syn2'),
                     input.id = "Input",
                     enrich.id = "Enrich")
# run ---- 
Enone <- enONE(Enone, 
               ruv.norm = TRUE, ruv.k = 3, 
               eval.pc.n = 5, eval.pam.k = 2:8,
               n.pos.eval = 200, n.neg.eval = 500)

## check performance ----
enScore <- getScore(Enone)
# perform PCA based on evaluation score, excluding BAT_SIM column (3) for no batch information provided, and SCORE column (9).
pca.eval <- prcomp(enScore[,-c(3, 9)], scale = TRUE)
PCA_Biplot(pca.eval, score = enScore$SCORE, pt.label = T)
# fig. 2E
ggsave('Norm_performance.pdf', width=8, height=6)

# select the top-ranked normalization
top.norm <- rownames(enScore[1,])
# [1] "TMM_RUVse_k3"

# perform top-ranked normalization on counts from dog (denoted as sample)
Enone <- UseNormalization(Enone, slot = 'sample', method = top.norm)
# get normalized counts from dog
top.norm.data <- Counts(Enone, slot = 'sample', method = top.norm)

# Find Enrichment ----
# NAD-RNAs were defined as fold change of normalized transcript counts ≥ 2, FDR < 0.05, 
# and log2-CPM > 1 in enrichment samples compared to those in input samples.
Enone <- FindEnrichment(Enone, slot = 'sample', method = top.norm)

# get filtered enrichment results 
res.sig.ls <- getEnrichment(Enone, slot='sample', filter=TRUE)
# further filtered by expression level
unlist(lapply(res.sig.ls, function(x) nrow(subset(x, logCPM>1))))
# D0_ctrl.Enrich_D0_ctrl.Input D14_ctrl.Enrich_D14_ctrl.Input 
# 2066                         1994                            
# D0_nmn.Enrich_D0_nmn.Input  D14_nmn.Enrich_D14_nmn.Input 
# 2545                        2184 

res.sig.ls <- lapply(res.sig.ls, function(x) subset(x, logCPM>1))

# dog nad-rna
# global dynamics of NAD-RNA by violin-box plot 
nad_df1 <- reduceRes(res.sig.ls, logfc.col = 'logFC')
nad_df1$Group <- gsub('\\..*', '', nad_df1$Group) # simplify group id
nad_df1$Group <- factor(nad_df1$Group, levels = unique(nad_df1$Group)) # order group id
BetweenStatPlot(nad_df1, x='Group', y='logFC', color='Group', step.increase = 0.6)
ggsave('results/StatVis/NAD-RNA_Global_bestNorm.pdf', width=10, height=6)

# save Enone data
save(Enone, file="Data/Enone.RData")

