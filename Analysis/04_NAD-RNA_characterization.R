# NAD-RNA characterization
library(enONE)
library(paintingr)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(webr)

##--Custom Function--##
# reduce un-filtered results list 
reduceResAll <- function(res.ls, logfc.col, levels=names(res.ls)) {
  df <- data.frame()
  for (id in names(res.ls)) {
    curr <- res.ls[[grep(id, names(res.ls), value=TRUE)]] 
    df1 <- curr %>% 
      dplyr::mutate(Group = factor(rep(id, nrow(curr)), levels = levels)) %>% 
      dplyr::select(GeneID, !!sym(logfc.col), logCPM, FDR, Group)
    df <- rbind(df, df1)
  }
  return(df)
}

# MA plot for NAD-RNA 
MADplot <- function(data, title=NULL) {
  ggdat <- data
  ggdat$sig <- ifelse(data$logFC >= 1 & data$FDR < 0.05 & data$logCPM >= 1, 
                      'NAD-RNA', 'None')
  ggdat <- ggdat[order(ggdat$sig, decreasing = TRUE),]
  
  x.pos <- max(ggdat$logCPM)*0.8
  y.pos <- max(ggdat$logFC)*0.9
  
  ggplot(ggdat, aes(logCPM, logFC, color=sig)) +
    geom_point(shape = 19) +
    annotate('text', x = x.pos, y = y.pos, 
             label = paste0('bold("NAD-RNA: ', sum(ggdat$sig == 'NAD-RNA'), '")'),
             color = '#046C9A', size = 3, parse = TRUE) +
    theme_classic() +
    theme(axis.text = element_text(color='black'),
          axis.ticks = element_line(color='black'),
          legend.position = 'none',
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 1, lty='dashed') +
    scale_color_manual(values = c("#046C9A", "#D5D5D3")) +
    labs(x='Log2 Average Abundance', y='Log2 Fold Change\n(Enrichment/Input)', 
         color='', title=title)
}

# wrapper of PieDonut for only showing donuts
donutChart <- function(data, title=NULL) {
  explode.donuts <- which(data$n <= 10)
  webr::PieDonut(data, 
           aes(Group, label, count=n),
           explode = 1, explodeDonut = T, r0 = 0.7, r1=0.7, 
           addPieLabel = F, showRatioPie = F, showRatioDonut = F, showRatioThreshold = 0.00001,
           selected = explode.donuts) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = paint_palette('Twilight',nrow(data),'continuous')) +
    labs(x='', y='', title = title)
}
##----##

# load NAD-RNA data from results of 03_enONE.R
load('Data/Enone.RData')
# get differential analysis results for all genes 
res.ls <- getEnrichment(Enone, slot='sample', filter=F)
# simplify group id
names(res.ls) <- c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')
nad_all <- reduceResAll(res.ls, logfc.col = 'logFC')

# get enriched genes
res.sig.ls <- getEnrichment(Enone, slot='sample', filter=T)
res.sig.ls <- lapply(res.sig.ls, subset, logCPM>1)
# simplify group id
names(res.sig.ls) <- c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')
nad_sig <- reduceRes(res.sig.ls, logfc.col = 'logFC')

# load gtf
## annotation retrieved from https://ftp.ensembl.org/pub/release-106/gff3/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf.gz
gtf.cf <- rtracklayer::import('Canis_lupus_familiaris.ROS_Cfam_1.0.106.gtf')
gene.cf <- as.data.frame(subset(gtf.cf, type == 'gene'))
anno.cols <- c('gene_id', 'gene_name', 'seqnames', 'gene_biotype') # only used these columns 

# merge NAD-RNA table with annotation
nad_sig <- nad_sig %>% left_join(gene.cf[, anno.cols], by=c('GeneID'='gene_id'))

# assign all chromosome start with JAA as unknown chromsome (Un) 
nad_sig$seqnames <- gsub('^JAA.*', 'Un', nad_sig$seqnames)
nad_sig$seqnames <- factor(nad_sig$seqnames, levels = c(as.character(1:38), 'X' ,'Un')) # order chr

# rename all *pseudogene as pseudogene
# rename IG_C_gene,IG_V_gene as IG_gene
# rename TR_V_gene,TR_C_gene as TR_gene
nad_sig$gene_biotype <- str_replace_all(nad_sig$gene_biotype, c('.*pseudogene'='pseudogene', '^IG.*'='IG_gene', '^TR.*'='TR_gene'))

# MA plot ----
# fig. 3A
# draw MA plot
madp.ls <- lapply(names(res.ls), function(i) {
  nad_all %>% 
    filter(Group == i) %>% 
    MADplot(title = i)
})
names(madp.ls) <- names(res.ls)
madps <- wrap_plots(madp.ls, nrow=1) & scale_y_continuous(limits=c(-15, 8), expand=expansion(c(0,0)))
ggsave('MAD_all.pdf', madps, width=24, height=5)

# gene type ----
# fig. 3B
df1 <- nad_sig %>% 
  dplyr::count(Group, gene_biotype) %>% 
  group_by(Group) 
# label each pie with gene biotype and numbers, e.g., protein-encoding (1951)
df1$label <- paste0(df1$gene_biotype, '\n', '(', df1$n, ')')
df1$Group <- as.character(df1$Group)
# draw donut charts
dp.ls <- lapply(names(res.ls), function(i) {
  df1 %>% 
    filter(Group == i) %>% 
    donutChart(title = i)
})
names(dp.ls) <- names(res.ls)
dps1 <- wrap_plots(dp.ls, nrow=1)
ggsave('GeneType_donuts.pdf', dps1, width=20, height=6)

# chromosome distribution ----
# calculate proportion of NAD-RNA by chromosome
df2 <- nad_sig %>% 
  dplyr::count(Group, seqnames) %>% 
  group_by(Group) %>% 
  mutate(pct = 100*n/sum(n)) 

# table for chromosome label 
label_data <- subset(df3, Group=='D0.Water')

# fig. 3C
# circular bar plot
ggplot(df2, aes(seqnames, pct)) +
  geom_hline(yintercept = c(0,10,20,30), color="#A8BAC4", size=0.3, linetype="solid") +
  geom_bar(aes(fill=Group),stat='identity') +
  geom_text(data=label_data, aes(x = seqnames, y = -1, label = seqnames),
            color="black", fontface="bold", size=3, inherit.aes=FALSE) +
  ylim(-25,30) +
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
  annotate(x = 40, y = 10, label = "10%", geom = "text", fontface="bold") +
  annotate(x = 40, y = 20, label = "20%", geom = "text", fontface="bold") +
  annotate(x = 40, y = 30, label = "30%", geom = "text", fontface="bold") +
  scale_fill_manual(values = c('#56565C','#A9A89F','#176AA6','#4F99B4')) +
  labs(fill='')
ggsave('ChrDistribution_circular_v2.pdf', width=10, height=8)


# gene length ~ logFC ----
# retrieve gene length from featrueCounts files
Translen <- read.table('Data/C1_counts.txt', sep='\t', header=T)

df3 <- left_join(nad_sig, Translen[,c('Geneid','Length')], by = c('GeneID'='Geneid')) %>% # merge with gene length info
  group_by(Group) %>% 
  mutate(quartile=factor(ntile(logFC, 5), levels = 1:5)) %>% # divide NAD-RNA by 5 declies 
  ungroup()

# fig. 3D
# draw boxplot
ggplot(df3, aes(quartile, Length/1000)) +
  geom_boxplot(aes(fill=quartile), color='white', outlier.color='black') +
  stat_summary(fun=mean, geom='line', aes(group=1), size=1, color='black') +
  stat_summary(fun=mean, geom='point', color='white') +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color='black')) +
  facet_wrap(~Group, nrow=2) +
  scale_fill_manual(values=paint_palette('Twilight',5,'continuous'), guide='none') +
  labs(x='Deciles', y='Gene Length (kb)')
ggsave('GeneLength_logFC.pdf', width=8, height=3)
