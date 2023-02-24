# dynamics of NAD-RNA
library(enONE)
library(tidyverse)
library(paintingr)
library(patchwork)
library(ComplexHeatmap)
library(biomaRt)

# Gene Trajectory plot
plotTrajectory <- function(data.plot, x, y, group, avg.span=2, title=NULL) {
  
  stat_df <- data.plot %>% 
    group_by(sample.id) %>% 
    mutate(med = median(!!sym(y)),
           avg = mean(!!sym(y))) %>% 
    ungroup() %>% 
    distinct(sample.id, .keep_all = T)
  
  gg1 <- ggplot(data.plot, aes_string(x=x, y=y)) +
    stat_smooth(geom="line", mapping = aes_string(group=group), 
                se=FALSE, span=2, color="#A9A9A0", alpha=0.3) +
    geom_smooth(data=stat_df, mapping=aes_(x = as.name(x), y = ~avg), 
                color="#1B6195", size=2,
                se=TRUE, method="loess", span=avg.span, fill="#5B8DB2") +
    annotate(geom="text", x = 25, y = -1.4, 
             label=paste0("n = ", length(unique(unlist(data.plot[, group]))))) +
    theme_classic() +
    theme(axis.text = element_text(color="black")) +
    labs(x="Age", y="Z-score", title=title)
  
  gg1 
}


# load data, including Enone, top.norm.data, and res.sig.ls (NAD-RNAs)
load("DATA/Enone.RData")

# In Data folder 
metadata <- read.csv("Data/metadata.csv")
# only use discovery cohort
metadata <- metadata[metadata$cohort.group == "Discovery",]
metadata.uni <- metadata[!duplicated(metadata$sample.id),]

# get annotation
hg_id <- grep("^ENSG", rownames(Enone), value=T)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_anno <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description", "entrezgene_id"),
                   values = hg_id,
                   filters = "ensembl_gene_id", 
                   mart = ensembl)

gene_anno$entrezgene_id[is.na(gene_anno$entrezgene_id)] <- ""
gene_anno <- gene_anno[!duplicated(gene_anno$ensembl_gene_id),]
rownames(gene_anno) <- gene_anno$ensembl_gene_id

# individual enrichment profiles ----
# log transformation of normalized counts
counts.norm.log <- log2(top.norm.data + 1)
enrich_idx <- CreateGroupMatrix(Enone$enrich) # first row is the index of enrichment samples and the second row is the index of input samples
ind.lfc <- counts.norm.log[, enrich_idx[1,]] - counts.norm.log[, enrich_idx[2,]]
colnames(ind.lfc) <- metadata.uni$sample.id

# get all NAD-RNA
nad_id <- unique(Reduce(union, lapply(res.sig.ls, function(x) x$GeneID)))

# keep NAD-RNA
## fold change >= 2
ind.lfc.keep.fc2 <- ind.lfc[nad_id_fc2,]
ind.scale.fc2 <- t(scale(t(ind.lfc.keep.fc2))) # scale

# expression profiles ----
expr.norm <- counts.norm.log[,enrich_idx[2,]]
colnames(expr.norm) <- meta.uni$sample.id
expr.scale <- t(scale(t(expr.norm))) # scale

# vis heatmap of NAD-RNA ----
# heatmap parameter
col_breaks <- c(-2,-1.5,0,1,3,4)
my_pal <- c("#053264","#569EC9","#FFFFFF","#F8B597","#C13438","#6C0321")
col_fun <- circlize::colorRamp2(col_breaks, my_pal)

# palette for age
age_col1 <- circlize::colorRamp2(c(20, 45, 70), c("#edd0ca","#aa688f","#2d1d3d"))
color_map <- age_col1(meta.uni$age)

# palette for gender
pal2 <- setNames(c("#000000", "#DCDDDD"), nm=c("M", "F"))

# individual profiles heatmap
column_ha1 <- HeatmapAnnotation(num = anno_barplot(colSums(ind.lfc.keep.fc2 >= 1), 
                                                   border = FALSE,
                                                   bar_width = 0.8,
                                                   gp = gpar(fill=color_map, lwd=NA),
                                                   add_numbers = TRUE,
                                                   numbers_rot = 90,
                                                   numbers_offset = unit(-0.8, "cm"),
                                                   numbers_gp = gpar(col="#FFFFFF"),
                                                   height = unit(3, "cm")),
                                age = metadata.uni$age,
                                gender = metadata.uni$sex,
                                col = list(age=color_map, gender=pal2))

# cutree with k=3
h.clust <- cutree(hclust(dist(ind.scale.fc2)), k = 3)

pdf(file = "results/Ind_NADRNA_profiles_age.pdf", width=16, height=10)
Heatmap(ind.scale.fc2,
        name = "Z-score", split=h.clust, clustering_method_columns = "ward.D2",
        show_column_names = F, show_row_names = F,
        col = col_fun, cluster_rows = T,
        column_names_rot = 45, column_names_centered = T,
        use_raster = F, top_annotation = column_ha1
)
dev.off()

# identification of age-associated trajectories ----
# nad-rna
h.clust.df <- data.frame(gene = names(h.clust),
                         cluster = h.clust,
                         row.names = names(h.clust))

nad_cluster <- as.data.frame(ind.scale.fc2) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = -gene, 
               names_to = "sample.id",
               values_to = "value") %>% 
  left_join(metadata.uni[,c("sample.id","age")], by="sample.id") %>% 
  left_join(h.clust.df, by="gene")

tp.lsi <- lapply(unique(nad_cluster$cluster), function(i) {
  plotTrajectory(subset(nad_cluster, cluster==i), x="age", y="value", 
                 group="gene", title=i, avg.span=1)
})
tpsi <- wrap_plots(tp.lsi, ncol=1)
ggsave(paste0("results/Trajectory_NAD_hc3.pdf"), tpsi, width=3, height=9)

# expression levels of genes from different NAD-RNA clusters
tp.ls.expr <- lapply(unique(nad_cluster$cluster), function(i) {
  plotTrajectory(subset(expr.long, gene %in% subset(nad_cluster, cluster==i)$gene), 
                 x="age", y="value", group="gene", avg.span=3)
})
tp.ls.expr <- wrap_plots(tp.ls.expr, nrow=1) + plot_annotation("Transcriptome")
ggsave(paste0("results/Trajectory_expr_NAD_hc3.pdf"), tp.ls.expr, width=9, height=3)

# compute correlation between age and nad modifition/gene expression
sample.age <- metadata.uni$age
# calculate spearman correlation coefficient between gene expression and ages
cor.expr.df <- apply(expr.norm, 1, function(x) {
  # Spearman"s correlation test
  test.i <- cor.test(x, sample.age, method="spearman")
  # result table
  c(cor=test.i$estimate, p=test.i$p.value)
})
cor.expr.df <- as.data.frame(t(cor.expr.df))
# adjust multiple testing by BH method
cor.expr.df$fdr <- p.adjust(cor.expr.df$p, method="fdr")

# calculate spearman correlation coefficient between nad-rna modification and ages
cor.nad.df <- apply(ind.lfc.keep.fc2, 1, function(x) {
  # Spearman"s correlation test
  test.i <- cor.test(x, sample.age, method="spearman")
  # result table
  c(cor=test.i$estimate, p=test.i$p.value)
})
cor.nad.df <- as.data.frame(t(cor.nad.df))
cor.nad.df$fdr <- p.adjust(cor.nad.df$p, method="fdr")

# save table
# expr
expr.cor <- cor.expr.df %>% 
  rownames_to_column("GeneID") %>% 
  left_join(gene_anno, by=c("GeneID"="ensembl_gene_id"))

# nad-rna
nad.cor <- cor.nad.df %>% 
  rownames_to_column("GeneID") %>% 
  left_join(gene_anno, by=c("GeneID"="ensembl_gene_id"))

write.csv(expr.cor, "results/AgeCor_Expr.csv", row.names=F, quote=F)
write.csv(nad.cor, "results/AgeCor_NAD.csv", row.names=F, quote=F)

# draw distribution of correlation by clusters
cor.df.combine <- data.frame(cor=c(cor.nad.df$cor.rho, cor.expr.df[nad_id,]$cor.rho),
                      type=rep(c("nad","expr"), each=length(nad_id)),
                      cluster=rep(h.clust.df[nad_id,]$cluster, times=2))

# perform KS test
ks.p <- sapply(1:3, function(i) {
  ks.test(cor~type, data=subset(cor.df.combine, cluster==i), exact=FALSE)$p.value
})

# compute the average correlation by cluster on transcriptome and epitranscriptome, respectively
mean_dat <- cor.df.combine %>% group_by(cluster, type) %>% summarise(m=mean(cor))

# draw distribution charts
cor.df.combine %>% 
  ggplot(aes(cor, fill=type)) +
  geom_density(alpha=0.6, color=NA) +
  geom_vline(aes(xintercept = m, color=type), data=mean_dat, lty="dashed") +
  facet_wrap(~cluster) + 
  theme_classic() +
  theme(axis.text = element_text(color="black"),
        strip.background = element_blank(),
        strip.text = element_text(face="bold")) +
  scale_fill_manual(values=rev(c("#C35743", "#191919"))) +
  scale_color_manual(values=rev(c("#C35743", "#191919"))) +
  scale_y_continuous(expand = expansion(0)) +
  labs(x="Spearman"s correlation with age", y="Density", color="", fill="")
ggsave("results/AgeCorDist.pdf", width=8, height=3)

# draw single gene trajectory ----
# interested genes
itg1 <- c("ENSG00000167004","ENSG00000151092","ENSG00000164305","ENSG00000064012",
          "ENSG00000151461","ENSG00000154146","ENSG00000105329")

# transfer into long format
# modification levels of interested genes
nad.long <- as.data.frame(ind.scale.fc2[itg1,]) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = -gene, 
               names_to = "sample.id",
               values_to = "value") %>% 
  left_join(meta.uni[,c("sample.id","age")], by="sample.id") 
nad.long$type <- rep("Epitranscriptome", nrow(nad.long))

# transfer into long format
# expression levels of interested genes
expr.long <- as.data.frame(expr.scale[itg1,]) %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = -gene, 
               names_to = "sample.id",
               values_to = "value") %>% 
  left_join(meta.uni[,c("sample.id","age")], by="sample.id") 
expr.long$type <- rep("Transcriptome", nrow(expr.long))

# combine
df.long.it <- rbind(nad.long, expr.long) %>% 
  left_join(gene_anno, by=c("gene"="ensembl_gene_id")) %>% 
  mutate(external_gene_name=factor(external_gene_name, levels=gene_anno[itg1,]$external_gene_name))

# draw trajectory by genes
df.long.it %>% 
  ggplot(aes(age, value, color=type)) +
  geom_smooth(aes(group=type), se=T) +
  theme_classic() +
  theme(axis.text = element_text(color="black"),
        strip.background = element_blank(),
        strip.text = element_text(face="bold")) +
  facet_wrap(~external_gene_name, nrow=1) +
  scale_color_manual(values=c("#C35743", "#191919")) +
  labs(x="Age (years)", y="Z-score", color="")
ggsave("results/Trajectory_selGene_combine.pdf", width=16, height=3)

# save data
save(nad_id, ind.lfc.keep.fc2, ind.scale.fc2, expr.norm, expr.scale, file="DATA/NADRNA_profiles.RData")