# Sequencing stat
# load required libraries
library(tidyverse)
library(paintingr)

# read in metadata
metadata <- read.csv("DATA/metadata.csv", header=TRUE, comment.char = "#")

# read in counts table
counts_df <- read.csv("DATA/Counts.csv", row.names = 1)

# prefix of Drosophila spike-in gene id
spikeIn <- "^FB"

## uniquely mapped read pairs of each library ----
align_stat <- read.csv("DATA/align_stat_all.csv")
# first column is the sample id, second column is the library ID
# and the third column is the number of uniquely mapped read pairs of each library
## supp. fig. 1A ----
## draw bar plot of sum of the alignments
mapbp1 <- align_stat %>% 
  ggplot(aes(ID, Uni/1e6)) +
  geom_col(aes(fill=batch), width=0.8) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top") +
  scale_fill_manual(values=paint_palette("Twilight")) +
  scale_y_continuous(expand = expansion(0)) +
  labs(x="Samples", y="Uniquely mapped read pairs (M)", fill="Batch")
ggsave("results/MappedReads.pdf", mapbp1, width=10, height=6)

# sequencing saturation ----
geneNums_df <- read.csv("DATA/geneNums_all.csv")
## supp. fig. 2B ----
## draw sequencing saturation charts
seqp1 <- geneNums_df %>% 
  ggplot(aes(x=ReadsM, y=Num)) +
  geom_line(aes(color=batch, group=Sample),size=0.8) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line("black"),
        strip.text = element_text(face="bold"),
        axis.text = element_text(color="black"),
        legend.position = "top") +
  facet_wrap(~factor(Assay, levels=c("Input","Enrich")), ncol=1) +
  scale_color_manual(values=paint_palette("Twilight")) +
  labs(x="Reads (M)", y = "Numbers of Mapped Genes", color="Batch")
ggsave("results/SeqSaturation.pdf", seqp1, width = 4, height = 6)

