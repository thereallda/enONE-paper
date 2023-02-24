# enONE for data normalization and NAD-RNA identification
library(enONE)
library(tidyverse)

# In Data folder 
metadata <- read.csv("Data/metadata.csv")
counts_df <- read.csv("Data/Counts.csv", row.names = 1)

# only use discovery cohort
metadata <- metadata[metadata$cohort.group == "Discovery",]
counts_df <- counts_df[, metadata$id]

# prefix of Drosophila spike-in gene id
spikeInPrefix <- "^FB"

## filtering ----
## filter low-expressed genes
keep <- edgeR::filterByExpr(counts_df, group = metadata$condition, min.count = 20)
counts_keep <- counts_df[keep,]

## filter rRNA and TEC of human
## annotation retrieved from https://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.94.gtf")
gene <- as.data.frame(subset(gtf, type == "gene"))
rrna.id <- subset(gene, gene_biotype == "rRNA")$gene_id
tec.id <- subset(gene, gene_biotype == "TEC")$gene_id
counts_keep <- counts_keep[!rownames(counts_keep) %in% c(rrna.id, tec.id), ]

# create Enone ----
Enone <- createEnone(counts_keep, 
                     bio.group = metadata$condition,
                     enrich.group = metadata$Assay,
                     batch.group = metadata$batch,
                     spike.in.prefix = spikeInPrefix,
                     synthetic.id = c("Syn1","Syn2"),
                     input.id = "Input",
                     enrich.id = "Enrich")

## outlier assessment
OutlierTest(Enone, return=FALSE)
## all samples pass the test

## run ---- 
Enone <- enONE(Enone, 
               ruv.norm = TRUE, ruv.k = 5, 
               eval.pc.n = 5, eval.pam.k = 2:12,
               return.norm = TRUE)

## check performance
enScore <- getScore(Enone)

# perform PCA based on evaluation score, excluding SCORE column (9).
pca.eval <- prcomp(enScore[,-c(9)], scale = TRUE)
PCA_Biplot(pca.eval, score = enScore$SCORE, pt.label = T)
# fig. 1c
ggsave("results/Norm_performance.pdf", width=8, height=6)

# select the top-ranked normalization
top.norm <- rownames(enScore[1,])
# [1] "DESeq_RUVg_k4"

# get normalized counts from human (denoted as sample)
top.norm.data <- Counts(Enone, slot = "sample", method = top.norm)

# PCA before and after normalization
p1 <- PCAplot(Counts(Enone, "sample", "Raw"), 
              color = Enone$batch,
              shape = Enone$enrich,
              palette = paint_palette("Twilight"),
              vst.norm = TRUE, title="Raw")

p2 <- PCAplot(log1p(top.norm.data), 
              color = Enone$batch,
              shape = Enone$enrich,
              palette = paint_palette("Twilight"),
              vst.norm = FALSE, title="Normalized")

# ggExtra::ggMarginal for adding margins
mp1 <- ggExtra::ggMarginal(p1, groupColour = TRUE, groupFill = TRUE)
mp2 <- ggExtra::ggMarginal(p2, groupColour = TRUE, groupFill = TRUE)
ggsave("results/PCA_sample_mar_raw.pdf", mp1, width=4,height=4)
ggsave("results/PCA_sample_mar_norm.pdf", mp2, width=4,height=4)

# Find Enrichment ----
# NAD-RNAs were defined as fold change of normalized transcript counts ≥ 2, FDR < 0.05, 
# and log2-CPM > 1 in enrichment samples compared to those in input samples.
Enone <- FindEnrichment(Enone, slot = "sample", method = top.norm)

# get filtered enrichment results 
res.sig.ls <- getEnrichment(Enone, slot="sample", filter=TRUE)
# further filtered by expression level
res.sig.ls <- lapply(res.sig.ls, function(x) subset(x, logCPM>=1))

# save data
save(Enone, res.sig.ls, top.norm.data, file="Data/Enone.RData")