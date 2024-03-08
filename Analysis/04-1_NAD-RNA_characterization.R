# NAD-RNA characterization
library(enONE)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(paintingr)

# load data, including Enone, top.norm.data, and res.sig.ls (NAD-RNAs)
load("Data/Enone.RData")

# get all genes 
res.ls <- getEnrichment(Enone, slot="sample", filter=FALSE)
# name of the elements
names(res.ls) <- c("Young", "Mid", "Old")

# name of the elements
names(res.sig.ls) <- c("Young", "Mid", "Old")
nad_sig <- reduceRes(res.sig.ls, logfc.col = "logFC")

# load gtf
## annotation retrieved from https://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.94.chr.gtf")
gene <- as.data.frame(subset(gtf, type == "gene"))
anno.cols <- c("gene_id", "gene_name", "seqnames", "gene_biotype")

# merge NAD-RNA table with annotation
nad_sig <- nad_sig %>% left_join(gene[, anno.cols], by=c("GeneID"="gene_id"))

# rename all *pseudogene as pseudogene
# consider IG_C_gene,IG_V_gene as IG_gene
# consider TR_V_gene,TR_C_gene as TR_gene
# consider bidirectional_promoter_lncRNA as lncRNA
nad_sig$gene_biotype <- str_replace_all(nad_sig$gene_biotype, 
                                        c(".*pseudogene"="pseudogene", 
                                          "^IG.*"="IG_gene", 
                                          "^TR.*"="TR_gene",
                                          "bidirectional_promoter_lncRNA"="lncRNA"))
# gene type ----
nad_sig$gene_biotype <- factor(nad_sig$gene_biotype,levels=unique(nad_sig$gene_biotype))

# compute percentage of each gene type
df1 <- nad_sig %>% 
  dplyr::count(gene_biotype, .drop=FALSE) %>% 
  mutate(pct=round(n/sum(n)*100,2),
         gene_biotype=forcats::fct_reorder(gene_biotype, pct, .desc = TRUE))

# barplot for genebiotype with y-axis break
bp1 <- ggplot(df1, aes(gene_biotype, pct)) +
  geom_bar(stat="identity", width = 0.8) +
  geom_text(aes(label=pct), data=subset(df1, pct<3), vjust=-1) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  scale_fill_manual(values = rev(paint_palette("Splash",3,"continuous"))) +
  scale_x_discrete(labels = gsub("_"," ",levels(df1$gene_biotype))) +
  labs(x="", y="Percentage (%)", fill="") +
  ggbreak::scale_y_break(c(20,59), scales = 0.3)
ggsave("results/GeneType_bar.pdf", bp1, width=8, height=6, onefile=FALSE)

# chromosome distribution ----
# compute proportion of chromosomes
df2 <- nad_sig %>% 
  dplyr::count(seqnames, .drop=TRUE) %>% 
  mutate(pct=round(n/sum(n)*100,2))

# circular bar plot
cbp2 <- ggplot(df2, aes(seqnames, pct)) +
  geom_hline(yintercept=seq(0,10,2.5), color="#A8BAC4", size=0.3, linetype="solid") +
  geom_bar(stat="identity", width=0.6, fill="#3381CD") +
  geom_text(aes(x = seqnames, y = -1, label = seqnames),
            color="black", fontface="bold", size=3, inherit.aes=FALSE) +
  ylim(-8,10) +
  coord_polar(start = 0) +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank()
  ) +
  # Annotate custom scale inside plot
  annotate(x = 22, y = 10, label = "10%", geom = "text", fontface="bold") +
  annotate(x = 22, y = 5, label = "5%", geom = "text", fontface="bold") 
ggsave("results/ChrDistr_circ.pdf", cbp2, width=8, height=6)

# gene length ~ logFC ----
# retrieve gene length from featrueCounts files
Translen <- read.table("DATA/B11_Counts.txt", sep="\t", header=T)

# group NAD-RNA into five declies based on enrichment 
df3 <- left_join(nad_sig, Translen[,c("Geneid","Length")], by = c("GeneID"="Geneid")) %>% 
  mutate(quartile=factor(ntile(logFC, 5), levels = 1:5)) 

bxp3 <- ggplot(df3, aes(quartile, Length/1000)) +
  geom_boxplot(aes(fill=quartile), color="white", outlier.color="black", width=0.7) +
  stat_summary(fun=mean, geom="line", aes(group=1), size=1, color="black") +
  stat_summary(fun=mean, geom="point", color="white") +
  theme_classic() +
  theme(axis.text = element_text(color="black")) +
  scale_fill_manual(values=paint_palette("Twilight",5,"continuous"), guide="none") +
  labs(x="Deciles", y="Gene Length (kb)")

# number of introns ~ logFC ----
Translen$nInt <- sapply(1:nrow(Translen), function(i) {
  length(str_split(Translen$Chr[i], ";", simplify = TRUE)) - 1
})

# group NAD-RNA into five declies based on enrichment 
df4 <- left_join(nad_sig, Translen[,c("Geneid","nInt")], by = c("GeneID"="Geneid")) %>% 
  mutate(quartile=factor(ntile(logFC, 5), levels = 1:5)) 

bxp4 <- ggplot(df4, aes(quartile, nInt)) +
  geom_violin(aes(fill=quartile), color=NA) +
  geom_boxplot(fill="black", outlier.color=NA, width=0.2) +
  stat_summary(fun=mean, geom="line", aes(group=1), size=1, color="black") +
  stat_summary(fun=mean, geom="point", color="white") +
  theme_classic() +
  theme(axis.text = element_text(color="black"),
        legend.position = "none") +
  scale_fill_manual(values=paint_palette("Twilight",5,"continuous")) +
  labs(x="Deciles", y="Number of Introns")

bxp3 + bxp4 + plot_layout(ncol=1)
ggsave("results/GeneLength_NumIntrons_logFC.pdf", width=2, height=5)
