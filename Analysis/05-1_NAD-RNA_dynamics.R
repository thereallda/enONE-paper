# dynamics of NAD-RNA
library(enONE)
library(tidyverse)
library(patchwork)
library(paintingr)
library(ComplexHeatmap)
library(rstatix)

##--Custom Function--##
# dual-barplot ----
dualBarplot <- function(enrich.obj, showCategory=10){
  eobj <- enrich.obj[,c('Description','Count','p.adjust')]
  eggtmp <- eobj %>% 
    head(.,n=showCategory) %>% 
    dplyr::mutate(p.adjust = -log10(p.adjust)) %>% 
    reshape2::melt()
  bp <- ggplot(eggtmp,aes(x=factor(Description, levels=rev(unique(Description))))) +
    geom_bar(aes(y=value, fill=variable), stat='identity', 
             position=position_dodge(width=0.9), width=0.8) +
    coord_flip() +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), lty=4, col="#C35743", lwd=0.6) +
    theme(axis.text = element_text(color='black'),
          axis.line = element_line(),
          axis.ticks = element_line(color='black'),
          legend.position = 'top',
          panel.grid = element_blank()) +
    scale_y_continuous(expand = expansion(add=c(0,0))) +
    scale_fill_manual(values = c("#13559F","#DF832C"), labels = c('Gene Counts', '-log10(p.adjust)')) +
    labs(x='',y='',fill='')
  return(bp)
}
# betweenstatplot only show box ----
BoxBetweenStatPlot <- function(data, x, y, color, palette = NULL,
                               test = c('wilcox.test', 't.test', 'none'),
                               comparisons = NULL,
                               step.increase=0.3) {
  stat.formula <- as.formula(paste(y, "~", x))
  
  test <- match.arg(test, choices = c('wilcox.test', 't.test', 'none'))
  if (test != 'none') {
    if (test == 'wilcox.test') {
      stat_dat <- data %>%
        rstatix::wilcox_test(stat.formula, comparisons = comparisons)
    }
    if (test == 't.test') {
      stat_dat <- data %>%
        rstatix::t_test(stat.formula, comparisons = comparisons)
    }
    stat_dat <- stat_dat %>%
      rstatix::p_format(p, digits = 2, leading.zero = FALSE,
               trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>%
      rstatix::add_xy_position(x = x, dodge=0.8, step.increase=step.increase)
  }
  
  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(as.factor(data[,x])),")")
  x.num <- length(unique(data[,color])) # number of x types
  if (is.null(palette)) palette <- paint_palette("Spring", x.num, 'continuous')
  
  p <- data %>%
    ggplot(aes_string(x, y, color = color)) +
    geom_boxplot(width = 0.6) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = x.labs) +
    labs(x='')
  
  if (exists('stat_dat')) {
    p <- p + ggpubr::stat_pvalue_manual(data = stat_dat, label = "p", tip.length = 0, size = 3)
  }
  
  return(p)
}
##----##

# load NAD-RNA data from results of 03_enONE.R
load('data/Enone.RData')

# get enriched genes
res.sig.ls <- getEnrichment(Enone, slot='sample', filter=T)
res.sig.ls <- lapply(res.sig.ls, subset, logCPM>1)
# simplify group id
names(res.sig.ls) <- c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')

# calculate scaled NAD-RNA modification matrix ----
# get sample raw coutns
counts_raw <- Counts(Enone, slot='sample', method='Raw')
prior.count <- mean(1e6*20/colSums(counts_raw)) # offset

# get normalized counts
counts_norm <- Counts(Enone, slot='sample', 'TMM_RUVse_k3')
enrich_idx <- CreateGroupMatrix(Enone$enrich)
counts_norm_log <- log2(counts_norm + prior.count)

# sample-wise log fold-change matrix
ind.fc <- counts_norm_log[, enrich_idx[1,]] - counts_norm_log[, enrich_idx[2,]]

# create sample name, e.g., D0_ctrl.Input.1
samples_name <- paste(Enone$condition, Enone$replicate, sep='.')
# reassign column names of matrix, e.g., D0_ctrl.1
colnames(ind.fc) <- unique(gsub('(Input.)|(Enrich.)', '', samples_name))

# get nad-rna id from all samples
nad.id <- purrr::reduce(lapply(res.sig.ls, function(x) x$GeneID), union)

# get nad-rna id from nmn samples
nad.id.nmn <- purrr::reduce(lapply(res.sig.ls[grep('NMN',names(res.ls))], function(x) x$GeneID), union)

# get nad-rna id from water control samples
nad.id.ctrl <- purrr::reduce(lapply(res.sig.ls[grep('Water',names(res.ls))], function(x) x$GeneID), union)

# define differential modified NAD-RNAs by group ----
# get logFC values of all genes in each group
lfc.avg <- purrr::reduce(lapply(getEnrichment(Enone, slot='sample', filter=F), function(x) x[,c('GeneID','logFC')]), left_join, by='GeneID')
rownames(lfc.avg ) <- lfc.avg $GeneID
colnames(lfc.avg ) <- c('GeneID','D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')
# comparing NAD level between D14 and D0 samples from NMN or control groups, respectively
lfc.avg  <- lfc.avg  %>% 
  dplyr::mutate(logfc.diff.ctrl=D14.Water-D0.Water,
                logfc.diff.nmn=D14.NMN-D0.NMN,
                species.ctrl=if_else(GeneID %in% nad.id.ctrl, 'NAD-RNA', 'Not'),
                species.nmn=if_else(GeneID %in% nad.id.nmn, 'NAD-RNA', 'Not'))

# differential modified NAD-RNA criteria
# control up-regulated NAD-RNAs
s1 <- subset(lfc.avg , logfc.diff.ctrl >= log2(1.5) & species.ctrl == 'NAD-RNA')$GeneID
# control down-regulated NAD-RNAs
s2 <- subset(lfc.avg , logfc.diff.ctrl <= log2(2/3) & species.ctrl == 'NAD-RNA')$GeneID
# NMN up-regulated NAD-RNAs
s3 <- subset(lfc.avg , logfc.diff.nmn >= log2(1.5) & species.nmn == 'NAD-RNA' & !GeneID %in% s1)$GeneID
# NMN down-regulated NAD-RNAs
s4 <- subset(lfc.avg , logfc.diff.nmn <= log2(2/3) & species.nmn == 'NAD-RNA' & !GeneID %in% s2)$GeneID

# heatmap palette ----
col_fun <- circlize::colorRamp2(breaks=c(-2,-1,0,1,2), 
                                colors=c("#000000","#121212","#1f1f1f","#D16D45","#FEB24C"))

# draw heatmap by group ----
# scale control logfc values
# shape of control scaled matrix: 2296 x 8
ind.fc.ctrl.scale <- t(scale(t(ind.fc[nad.id.ctrl, 1:8])))
# trimming values for heatmap visualization
ind.fc.ctrl.scale[ind.fc.ctrl.scale > 2] <- 2
ind.fc.ctrl.scale[ind.fc.ctrl.scale < -2] <- -2

# scale nmn logfc values
# shape of nmn scaled matrix: 2705 x 10
ind.fc.nmn.scale <- t(scale(t(ind.fc[nad.id.nmn, 9:18])))
# trimming values for heatmap visualization
ind.fc.nmn.scale[ind.fc.nmn.scale > 2] <- 2
ind.fc.nmn.scale[ind.fc.nmn.scale < -2] <- -2

# manually ordering differential NAD-RNA
# order differential NAD-RNAs by decreased NAD-RNA and increased NAD-RNA
nad.diff.id.ctrl <- c(s2, s1) 
nad.diff.id.nmn <- c(s4, s3)
# order non-differential NAD-RNAs by average hierachical clustering based on euclidean distances
nad.nondiff.id.ctrl <- nad.id.ctrl[!nad.id.ctrl %in% nad.diff.id.ctrl]
clust.ctrl <- hclust(dist(ind.fc.ctrl.scale[nad.nondiff.id.ctrl,]))$order
nad.nondiff.id.nmn <- nad.id.nmn[!nad.id.nmn %in% nad.diff.id.nmn]
clust.nmn <- hclust(dist(ind.fc.nmn.scale[nad.nondiff.id.nmn,]))$order

# define row order
row_ord_ctrl <- c(nad.diff.id.ctrl, nad.nondiff.id.ctrl[clust.ctrl])
# nmn row order
row_ord_nmn <- c(nad.diff.id.nmn, nad.nondiff.id.nmn[clust.nmn])

# row annotation
# for control heatmap
row_anno_ctrl <- rowAnnotation(set=c(rep("D14/D0_decreased",length(s4)), 
                                     rep("D14/D0_increased",length(s1)), 
                                     rep("Shared",length(nad.nondiff.id.ctrl))),
                               col=list(set=c("D14/D0_decreased"="#176AA6",
                                              "D14/D0_increased"="#C76041",
                                              "Shared"="#211F1E")))
# for nmn heatmap
row_anno_nmn <- rowAnnotation(set=c(rep("D14/D0_decreased",length(s3)), 
                                    rep("D14/D0_increased",length(s2)), 
                                    rep("Shared",length(nad.nondiff.id.nmn))),
                              col=list(set=c("NMN-down"="#176AA6",
                                             "D14/D0_increased"="#C76041",
                                             "Shared"="#211F1E")))
# column annotation
col_anno_ctrl <- HeatmapAnnotation(day = rep(c('D0','D14'), each=4),
                              col = list(day = c('D0'='#A9A9A0','D14'='#2D9E9E')))
col_anno_nmn <- HeatmapAnnotation(day = rep(c('D0','D14'), each=5),
                                   col = list(day = c('D0'='#A9A9A0','D14'='#2D9E9E')))
# ctrl heatmap
colnames(ind.fc.ctrl.scale) <- gsub('D.*\\.', '', colnames(ind.fc.ctrl.scale))
hm.ctrl <- Heatmap(ind.fc.ctrl.scale[row_ord_ctrl,],
        name='z-score', 
        col = col_fun,
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        use_raster = FALSE,
        row_split = c(rep("1",length(nad.diff.id.ctrl)), rep("2",length(nad.nondiff.id.ctrl))),
        row_title = ' ',
        column_title = 'Water',
        heatmap_legend_param = list(at = c(-2,-1,0,1,2)),
        top_annotation = col_anno_ctrl,
        left_annotation = row_anno_ctrl
)
# nmn heatmap
colnames(ind.fc.nmn.scale) <- gsub('D.*\\.', '', colnames(ind.fc.nmn.scale))
hm.nmn <- Heatmap(ind.fc.nmn.scale[row_ord_nmn,],
                   name='z-score', 
                   col = col_fun,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE,
                   use_raster = FALSE,
                   row_split = c(rep("1",length(nad.diff.id.nmn)), rep("2",length(nad.nondiff.id.nmn))),
                   row_title = ' ',
                   column_title = 'NMN',
                   heatmap_legend_param = list(at = c(-2,-1,0,1,2)),
                   top_annotation = col_anno_ctrl,
                   left_annotation = row_anno_nmn
)

# draw heatmap for control
pdf('results/StatVis/NAD_Global_heatmap_ctrl_v2.pdf', width=3, height=10)
draw(hm.ctrl)
dev.off()
# draw heatmap for NMN
pdf('results/StatVis/NAD_Global_heatmap_nmn_v2.pdf', width=3.15, height=10)
draw(hm.nmn)
dev.off()

# functional enrichment ----
library(clusterProfiler)
library(org.Cf.eg.db)
library(biomaRt)

# get all expressed gene id
all_dogid <- rownames(Counts(Enone, slot='sample', 'Raw'))

# gene annotation from ensembl
ensembl <- useMart("ensembl", dataset = "clfamiliaris_gene_ensembl")
gene_id <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description', 'entrezgene_id'),
                 values = all_dogid,
                 filters = 'ensembl_gene_id', 
                 mart = ensembl)
gene_id <- gene_id[!is.na(gene_id$entrezgene_id),]

# convert ensembl id to entrez id
convrt.ls <- lapply(res.ls, function(x) {
  x <- x %>% inner_join(gene_id, by=c('GeneID'='ensembl_gene_id'))
  return(x)
})
names(convrt.ls) <- c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')

# perform functional enrichment analysis based on NAD-RNAs from Day 0 dogs
ego1 <- enrichGO(union(convrt.ls$D0.ctrl$entrezgene_id, convrt.ls$D0.nmn$entrezgene_id),
                 keyType = "ENTREZID",
                 OrgDb = org.Cf.eg.db, ont = "BP", 
                 universe = as.character(gene_id$entrezgene_id), 
                 readable = TRUE)

# save table
write.csv(ego1@result, file='D0_GOBP.csv', row.names=F)

# simplify GO terms
ego1.sim <- simplify(ego1, cutoff=0.7)

# fig. 4B
# draw dual barplot
dualBarplot(ego1.sim[-c(2,5,12),])
ggsave('results/fe/D0_GOBP.pdf', width = 6, height = 3)

# comparing modification levels of NAD-RNAs from selected GO terms ----
# get differential analysis results for all genes 
res.all.ls <- getEnrichment(Enone, "sample")
# simplify names
names(res.all.ls) <- c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')

# reduce list into table 
nad_df <- reduceRes(res.all.ls, logfc.col = 'logFC')

# merge with ensembl annotation
nad_df <- left_join(nad_df, gene_id, by=c('GeneID'='ensembl_gene_id'))

# investigate following GO terms:
## GO:0006518 peptide metabolic process
## GO:0006412 translation
## GO:0007005 mitochondrion organization
## GO:0008380 RNA splicing
## GO:0022613 ribonucleoprotein complex biogenesis
## GO:0006955 immune response
## GO:0070661 leukocyte proliferation
## GO:0009060 aerobic respiration
## GO:0046034 ATP metabolic process
## GO:0006119 oxidative phosphorylation
selectGO <- data.frame(go=c('GO:0006518','GO:0006412','GO:0007005','GO:0008380','GO:0022613',
                            'GO:0006955','GO:0070661','GO:0009060','GO:0046034','GO:0006119'),
                       term=c('Peptide metabolism','Translation','Mitochondrion organization',
                              'RNA splicing','RNA-protein complex biogenesis',
                              'Immune response','Leukocyte proliferation',
                              'Aerobic respiration','ATP metabolism', 'Oxidative phosphorylation'))

## modification level boxplot by GO ----
# boxplot of NAD-RNA modification levels by group for each GO terms
# store in list
bp.ls1 <- list()
for (i in 1:nrow(selectGO)) {
  gogeneid1 <- unique(get(selectGO$go[i], org.Cf.egGO2ALLEGS))
  nad_sub_nmn <- subset(nad_df, entrezgene_id %in% gogeneid1 & GeneID %in% nad.id.nmn & Group %in% c('D0_nmn','D14_nmn'))
  nad_sub_ctrl <- subset(nad_df, entrezgene_id %in% gogeneid1 & GeneID %in% nad.id.ctrl & Group %in% c('D0_ctrl','D14_ctrl'))
  nad_sub <- rbind(nad_sub_ctrl,nad_sub_nmn)
  bp.ls1[[i]] <- BoxBetweenStatPlot(nad_sub, x='Group', y='logFC', color='Group', 
                            palette = c('#56565C','#A9A89F','#176AA6','#4F99B4'),
                            comparisons = list(c('D0_ctrl','D14_ctrl'), c('D0_nmn','D14_nmn')),
                            step.increase = 0.2) + ggtitle(selectGO$term[i])
}
# combine list of plots into one
bps1 <- wrap_plots(bp.ls1, nrow=2) & 
  scale_y_continuous(limits = c(-0.1,5.5)) & 
  theme(plot.title = element_text(face='bold', hjust=0.5, color="#C35743"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# remove redundant y-axis ticks, text and title 
for (i in c(2:5,7:10)) {
  bps1[[i]] <- bps1[[i]] + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())
}

# add x-axis text and ticks at bottom
for (i in 6:10) {
  bps1[[i]] <- bps1[[i]] + 
    theme(axis.ticks.x = element_line(color='black'),
          axis.text.x = element_text(color='black')) +
    scale_x_discrete(labels=c('D0','D14','D0','D14'))
}

# fig. 4C
# add legend
bps1 <- bps1 + theme(legend.position = 'right') + labs(color='')
ggsave('NAD_level_GO.pdf', bps1, width=8, height=4.5)

# expression levels of selected GO terms ----
# get normalized counts from dog
expr.mat <- as.data.frame(t(counts_norm[,enrich_idx[2,]]))
# average expression by group
expr.mat$Group <- Enone$condition[grep('Input', Enone$condition)]
avg.expr.mat <- aggregate(expr.mat[,-ncol(expr.mat)], list(expr.mat[,'Group']), FUN=mean)
avg.expr.mat$Group.1 <- NULL
avg.expr.mat <- as.data.frame(t(avg.expr.mat))
colnames(avg.expr.mat) <- c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN')
# log transformation with piror.count offset
avg.expr.mat <- log2(avg.expr.mat + prior.count)
# transform to long format
expr_df <- avg.expr.mat %>% 
  rownames_to_column("GeneID") %>% 
  pivot_longer(cols = -GeneID,
               names_to = 'Group',
               values_to = 'logexpr') %>% 
  left_join(gene_id, by=c('GeneID'='ensembl_gene_id'))
expr_df$Group <- factor(expr_df$Group, levels=c('D0.Water', 'D14.Water', 'D0.NMN', 'D14.NMN'))

## expression boxplot by GO ----
bp.ls2 <- list()
for (i in 1:nrow(selectGO)) {
  gogeneid1 <- unique(get(selectGO$go[i], org.Cf.egGO2ALLEGS))
  expr_sub_nmn <- subset(expr_df, entrezgene_id %in% gogeneid1 & GeneID %in% nad.id.nmn & Group %in% c('D0.nmn','D14.nmn'))
  expr_sub_ctrl <- subset(expr_df, entrezgene_id %in% gogeneid1 & GeneID %in% nad.id.ctrl & Group %in% c('D0.ctrl','D14.ctrl'))
  expr_sub <- as.data.frame(rbind(expr_sub_ctrl, expr_sub_nmn))
  bp.ls2[[i]] <- BoxBetweenStatPlot(expr_sub, x='Group', y='logexpr', color='Group', 
                                    palette = c('#56565C','#A9A89F','#176AA6','#4F99B4'),
                                    comparisons = list(c('D0.ctrl','D14.ctrl'), c('D0.nmn','D14.nmn')),
                                    step.increase = 0.2) + ggtitle(selectGO$term[i])
}

bps2 <- wrap_plots(bp.ls2, nrow=2) &
  scale_y_continuous(limits = c(0,15)) & 
  theme(plot.title = element_text(face='bold', hjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )
# remove redundant y-axis ticks, text and title 
for (i in c(2:5,7:10)) {
  bps2[[i]] <- bps2[[i]] + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())
}
# add x-axis text and ticks at bottom
for (i in 6:10) {
  bps2[[i]] <- bps2[[i]] + 
    theme(axis.ticks.x = element_line(color='black'),
          axis.text.x = element_text(color='black', angle=45, hjust=0.5, vjust=0.5)) +
    scale_x_discrete(labels=c('D0','D14','D0','D14'))
}

# supp. fig. 4B
# add legend
bps2 <- bps2 + theme(legend.position = 'right') + labs(color='')
ggsave('Expr_level_GO.pdf', bps2, width=10, height=6)

# save data for downstream analysiss
save(lfc.avg, ind.fc, nad_df, avg.expr.mat, expr_df, s1, s2, s3, s4, gene_id, file='Data/StatDat.RData')
