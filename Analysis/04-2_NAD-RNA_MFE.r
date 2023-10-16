# Analysis of MFE
library(enONE)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(paintingr)

# load data, including Enone, top.norm.data, and res.sig.ls (NAD-RNAs)
load("Data/Enone.RData")

nad_id <- unique(Reduce(union, lapply(res.sig.ls, function(x) x$GeneID)))

# load gtf
## annotation retrieved from https://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.94.chr.gtf")

# get coordination of 5' UTR and export as bed format ----
gtf_5utr <- as.data.frame(subset(gtf, type == "five_prime_utr"))
bed_5utr <- gtf_5utr[,c("seqnames","start","end","gene_id","score","strand")]

# for NAD-RNA
bed_5utr <- subset(bed_5utr, gene_id %in% nad_id)

# set region as the 5'utr start site + 100 bp downstream
bed_5utr$start <- ifelse(bed_5utr$strand == "+", bed_5utr$start, bed_5utr$end-99)
bed_5utr$end <- ifelse(bed_5utr$strand == "+", bed_5utr$start+99, bed_5utr$end)

# add prefix "chr"
bed_5utr$seqnames <- paste0("chr",bed_5utr$seqnames)
bed_5utr <- bed_5utr[!duplicated(bed_5utr$start),]
write.table(bed_5utr, file = "results/NADRNA_5UTR.bed", sep = "\t", col.names = F, quote = F, row.names = F, na=".")

# For non-NAD-RNA (n=600)
allgene <- rownames(Enone)
allgene <- grep("ENSG",allgene,value=T)

# randomly select 1000 times for permutation
# random 5utr
non_nad <- allgene[!allgene %in% nad_id]
non_nad <- non_nad[non_nad %in% unique(gtf_5utr$gene_id)]

bed_5utr_non_df <- data.frame()
for (i in 1:1000) {
  set.seed(i)
  non_nad1 <- sample(non_nad, size=600)
  
  bed_5utr_non <- gtf_5utr[,c("seqnames","start","end","gene_id","score","strand")]
  bed_5utr_non <- subset(bed_5utr_non, gene_id %in% non_nad1)
  
  # set region as the 5'utr start site + 100 bp downstream
  bed_5utr_non$start <- ifelse(bed_5utr_non$strand == "+", bed_5utr_non$start, bed_5utr_non$end-99)
  bed_5utr_non$end <- ifelse(bed_5utr_non$strand == "+", bed_5utr_non$start+99, bed_5utr_non$end)
  
  # add prefix "chr"
  bed_5utr_non$seqnames <- paste0("chr",bed_5utr_non$seqnames)
  bed_5utr_non <- bed_5utr_non[!duplicated(bed_5utr_non$start),]
  bed_5utr_non$seed <- rep(i, nrow(bed_5utr_non))
  bed_5utr_non_df <- rbind(bed_5utr_non_df, bed_5utr_non)
}
write.table(bed_5utr_non_df, file="results/NonNAD_futr_random1000.bed",sep = "\t", col.names = F, quote = F, row.names = F, na=".")

##----##
# RUN 04-2_RNAFold.sh in shell #
# which return MFE based on the sequence of 5'UTR#
##----##

# Visualization ----
# read MFE from RNAfold
# for NAD-RNA
nad_mfe <- read.csv("results/NADRNA_MFE.csv")

nad_mfe <- nad_mfe %>% 
  group_by(ID) %>% 
  summarise(mfe=mean(MEF))

# NAD-RNA into 5 declies based on enrichment
nad_sig <- reduceRes(res.sig.ls, logfc.col = "logFC")
nad_df <- nad_sig %>% 
  group_by(GeneID) %>% 
  summarise(logfc=log2(mean(2^logFC))) %>% 
  mutate(quartile=factor(ntile(logfc, 5), levels = 1:5)) 

df2 <- nad_mfe %>% 
  left_join(nad_df, by=c("ID"="GeneID"))

# for random genes
nonnad_mfe_rad <- read.csv("results/NonNAD_MFE_random1000.csv", skip=1)
colnames(nonnad_mfe_rad) <- c("GeneID", "mfe", "coord", "seed")
nonnad_mfe_rad$mfe <- as.numeric(nonnad_mfe_rad$mfe)
nonnad_mfe_rad <- na.omit(nonnad_mfe_rad)

nad_mfe_rad <- nonnad_mfe_rad %>% 
  group_by(seed, GeneID) %>% 
  summarise(mfe=mean(mfe))
nad_mfe_rad_avg <- nad_mfe_rad %>% 
  group_by(seed) %>% 
  summarise(avg_mfe=mean(mfe))
table(nad_mfe_rad_avg$avg_mef < mean(nad_mfe$mfe))

hp1 <- ggplot(nad_mfe_rad_avg, aes(avg_mfe)) +
  geom_histogram(fill="#bfbfbf") +
  geom_vline(xintercept = mean(nad_mfe$mfe), lty="dashed", color="#C35743") +
  annotate("text", x = -30.5, y = 100, color = "#C35743",
           label = "Average MFE of \ngenes with NAD-caps") +
  theme_classic() +
  theme(axis.text = element_text(color="black")) +
  labs(x="Average MFE (kcal/mol)", y="Frequency", subtitle = "P < 2e-16")

vp2 <- ggplot(df2, aes(quartile, mfe)) +
  geom_violin(aes(fill=quartile), width=0.7, color=NA) +
  stat_summary(fun = mean,
               fun.max = function(x) mean(x)+sd(x),
               fun.min = function(x) mean(x)-sd(x),
               geom = "pointrange", color = "white") +
  stat_summary(fun=mean, geom="line", aes(group=1), size=1, color="black") +
  theme_classic() +
  theme(axis.text = element_text(color="black")) +
  scale_fill_manual(values=paint_palette("Twilight",5,"continuous"), guide="none") +
  labs(x="Deciles", y="MFE", subtitle = "ANOVA, P = 1.1e-06")

hp1 + vp2


