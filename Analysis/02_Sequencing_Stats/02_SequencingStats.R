# Sequencing stat
# load required libraries
library(tidyverse)
library(paintingr)

# read in metadata
metadata <- read.csv("metadata.csv", header=TRUE, comment.char = "#")

# read in counts table
counts_df <- read.csv("Counts.csv", row.names = 1)

# prefix of Drosophila spike-in gene id
spikeIn <- "^FB"

# create sample name, e.g., D0_ctrl.Input.1
samples_name <- paste(metadata$condition, metadata$replicate, sep='.')

## uniquely mapped read pairs of each library ----
align_stat <- read.table('alignment_stat.tsv', sep='\t', header = TRUE)[,c(1,2)]
# first column is the library id and second column is the number of uniquely 
# mapped read pairs of each library
colnames(align_stat) <- c('ID', 'Uni')

# combine align_star and metadata
align_stat <- align_stat %>% 
  right_join(metadata, by=c('ID'='id')) %>% 
  mutate(group=str_extract(condition, '(ctrl)|(nmn)'),
         ID=factor(ID, levels=str_sort(ID, numeric = TRUE))) # sort id by numeric order
        
## supp. fig. 2A ----
## draw bar plot of sum of the alignments
align_stat %>% 
  ggplot(aes(ID, Uni/1e6)) +
  geom_col(aes(fill=group), width=0.8) +
  geom_text(aes(label=round(Uni/1e6,1)), color='white', nudge_y = -2) +
  coord_flip() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    axis.ticks.length = unit(0, "mm"),
    axis.line.y.left = element_line(color = "black"),
    axis.text = element_text(color = 'black'),
    legend.position = 'top'
  ) +
  scale_x_discrete(labels = str_remove(samples_name,'(ctrl).|(nmn).')) +
  scale_fill_manual(values=c('#56565C','#528EA1'), labels=c('Water','NMN')) +
  labs(x='', y='Uniquely mapped read pairs (M)', fill='')
ggsave('MappedReads.pdf', width=8, height=8)

# sequencing saturation ----
geneNums <- read.table('geneNum.txt', sep='\t',header=TRUE)
# first column is the subsampling proportion of library; 
# second column is the number of well-detected genes (counts>10) in corresponding proportion; 
# thrid column is the id of library
colnames(geneNums) <- c('Per','Num','Sample')

# get sum of the read counts derived from dog
libSize <- data.frame(libsize=colSums(counts_df[grep(spikeIn, rownames(counts_df), invert = TRUE),]),
                      Sample=colnames(counts_df))

# add library size 
geneNums_df <- geneNums %>% 
  # add zero coordinate
  add_row(Per=rep(0,length(unique(geneNums$Sample))), 
          Num=rep(0,length(unique(geneNums$Sample))), 
          Sample=unique(geneNums$Sample)) %>% 
  arrange(Sample, Per) %>% 
  left_join(libSize, by="Sample") %>% 
  mutate(ReadsM=round(Per*libsize/1e6, 3),
         Sample = factor(Sample, levels = unique(str_sort(Sample, numeric = TRUE)))) 

## supp. fig. 2B ----
# create labels with only day and replicate, e.g., D0.1
samples_label <- setNames(str_remove_all(samples_name,'(ctrl).|(nmn).|(Input).|(Enrich).'), nm=metadata$id)

## draw sequencing saturation charts
seqp1 <- geneNums_df %>% 
  ggplot(aes(x=ReadsM, y=Num)) +
  geom_line(aes(color=group, linetype=assay),size=0.8) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line('black'),
        strip.text = element_text(face='bold'),
        axis.text = element_text(color='black'),
        legend.position = 'top') +
  facet_wrap(~Sample, ncol=4, labeller = as_labeller(samples_label)) +
  scale_color_manual(values=c('#A9A9A0','#4F99B4'), labels=c('Control','NMN')) +
  labs(x='Reads (M)', y = 'Numbers of Mapped Genes', color='', linetype='')
ggsave('SeqSaturation.pdf', seqp1, width = 8, height = 6)

# fold change of synthetic RNA counts ----
# Statistics summary (mean and +/- sd)
.mean_sd <- function(x) {
  m <- mean(x)
  ymin <- m - stats::sd(x)
  ymax <- m + stats::sd(x)
  return(c(y=m, ymin=ymin, ymax=ymax))
}
##--##

# Syn1: synthetic spike-in RNA with 5% NAD-caps, which should be enriched
# Syn2: synthetic spike-in RNA with 100% m7G-caps, which should not be enriched
syn_counts <- counts_df[c('Syn1','Syn2'),]

# calculate fold change of synthetic RNA counts (raw) by enrichment/input
in.idx <- grep('Input',metadata$condition) # get column index of input
en.idx <- grep('Enrich',metadata$condition) # get column index of enrichment
syn_lfc <- log2(syn_counts[,en.idx]) - log2(syn_counts[,in.idx])

# transform to long format for ggplot
syn_df <- as.data.frame(syn_lfc) %>%
  rownames_to_column("syn_id") %>% 
  pivot_longer(cols = -syn_id,
               names_to = 'id',
               values_to = 'logfc') %>% 
  left_join(as.data.frame(metadata[,c('id', 'condition')]), by = 'id')
## supp. fig. 2C ----
## draw bar plot for fold change of synthetic RNA raw counts
syn_df %>% 
  ggplot(aes(syn_id, logfc)) +
  geom_bar(stat='summary', fun='mean', width=0.2, fill='grey70') +
  stat_summary(fun.data = .mean_sd, geom = 'errorbar', width=0.12, size=1) +
  geom_hline(yintercept = 0, color='black') +
  geom_hline(yintercept = 1, color='#C30D23', lty='dashed') +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color='black'),
        axis.text.x = element_text(angle=45, color='#0F3481', hjust=0.5, vjust=0.5)) +
  scale_x_discrete(label = c('Spike-in 2\n(5% NAD-RNA)', 
                             'Spike-in 3\n(100% m7G-RNA)')) +
  labs(x='', y='Log2 Fold Change')
ggsave('Raw_spikein_lfc.pdf', width=4, height=3)  
