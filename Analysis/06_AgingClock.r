# Age prediction model 
# features from 61 samples
# training with 35 samples
# validating with 26 samples
library(enONE)
library(tidyverse)
library(patchwork)
library(caret)
library(glmnet)
library(patchwork)
library(foreach)
library(parallel)
library(biomaRt)
library(randomForest)
source("helper.R")

# load data, including Enone, top.norm.data, and res.sig.ls (NAD-RNAs)
load("Data/Enone.RData")
# load data, including nad_id, ind.lfc.keep.fc2, ind.scale.fc2, expr.norm, expr.scale, cor.nad.df, cor.expr.df
load("Data/NADRNA_profiles.RData")

# In Data folder 
metadata <- read.csv("Data/metadata.csv")

# only use discovery cohort
meta.disc <- metadata[metadata$cohort.group == "Discovery",]
meta.uni <- meta.disc[!duplicated(meta.disc$sample.id),]
meta.uni$age_group <- gsub("\\..*", "", meta.uni$condition)

# modeling ----
## model based on epitranscriptome ----
# use age-associated NAD-RNAs
nad_cor <- subset(cor.nad.df, abs(cor.rho) > 0.3 & p < 0.05)$GeneID

# elastic net ----
# prepare data
model_df1 <- data.frame(t(ind.lfc.keep.fc2[nad_cor,meta.uni$sample.id]), age=meta.uni$age)

# fit model
model_nad <- glmNetFit(model_df1, y="age", seed=655, train.prop = 0.7, 
                       model = "glmnet", split = T, 
                       sample.group = meta.uni$age_group)

# get selected variables, excluding intercept (the first element)
sel_var_nad <- coef(model_nad$fit)[which(abs(as.vector(coef(model_nad$fit))) > 0),]

# LOOCV prediction 
nad_loo <- LOOP(model_df1, y="age", model = "glmnet", 
                alpha = model_nad$bestTune$alpha, 
                lambda = model_nad$bestTune$lambda)
# vis with scatter plot                 
loop_net1 <- modelScatter(prediction = nad_loo$prediction, evalutaion = nad_loo$evaluation,
                          title = paste0("Epitranscriptome (n = ", length(sel_var_nad)-1,")"))

## model based on transcriptome ----
# prepare data
# use age-associated genes
expr_cor <- subset(cor.expr.df, abs(cor.rho) > 0.5 & p < 0.05)$GeneID
model_df2 <- data.frame(t(expr.norm[expr_cor, meta.uni$sample.id]), age=meta.uni$age)

# fit model
model_expr <- glmNetFit(model_df2, y="age", seed=275, train.prop = 0.7, 
                        split = T, model = "glmnet",
                        sample.group = meta.uni$age_group)


# get selected variables, excluding intercept (the first element)
sel_var_expr <- coef(model_expr$fit)[which(abs(as.vector(coef(model.expr$fit))) > 0),]

# LOOCV prediction 
expr_loo <- LOOP(model_df2, y="age", model = "glmnet", 
                 alpha = model_expr$bestTune$alpha, 
                 lambda = model_expr$bestTune$lambda)
# vis with scatter                 
loop_net2 <- modelScatter(prediction = expr_loo$prediction, evalutaion = expr_loo$evaluation,
                          title=paste0("Transcriptome (n = ", length(sel_var_expr)-1,")"))


## model based on epitranscriptome and transcriptome ----
expr_sel <- t(expr.norm[names(sel_var_expr)[-1], meta.uni$sample.id])
colnames(expr_sel) <- paste0(colnames(expr_sel), ".expr")
nad_sel <- t(ind.lfc.keep.fc2[names(sel_var_nad)[-1], meta.uni$sample.id])
colnames(nad_sel) <- paste0(colnames(nad_sel), ".nad")


# combined model 
model_df3 <- data.frame(nad_sel, expr_sel, age=meta.uni$age)
model_combined <- glmNetFit(model_df3, y="age", seed=14, 
                            train.prop = 0.7, split = T, model = "glmnet",
                            sample.group = meta.uni$age_group)
sel_var_combined <- coef(model_combined$fit)[which(abs(as.vector(coef(model_combined$fit))) > 0),]
nad_var_num <- length(grep("nad",names(sel_var_combined[-1])))
expr_var_num <- length(sel_var_combined)-nad_var_num-1

combined_loo_el <- LOOP(model_df3, y="age", model = "glmnet", 
                        alpha = model_combined$bestTune$alpha, 
                        lambda = model_combined$bestTune$lambda)
loop_net3 <- modelScatter(prediction = combined_loo_el$prediction, 
                          evalutaion = combined_loo_el$evaluation,
                          title=paste0("Combined (n = ",expr_var_num,"+",nad_var_num,")"))

# Validation ----
meta.val <- subset(meta, cohort.group == "Validation")
meta.uni.val <- meta.val[!duplicated(meta.val$sample.id), ]
meta.uni.val$age_group <- gsub("\\..*","",meta.uni.val$condition)

# nad levels of validation cohort
ind.lfc.val <- ind.lfc.keep.fc2[,meta.uni.val$sample.id]

# expression levels of validation cohort
expr.val <- expr.norm[,meta.uni.val$sample.id]

## validation ----
# combined model: 31 = 20 + 11
valid.m <- validPred(nad.mat = ind.lfc.val, expr.mat = expr.val, 
                     model = model_combined$fit, obs = metadata.uni.val$age)

p7 <- modelScatter(valid.m$prediction, valid.m$evaluation, 
                   title="Validation (N=26)\nCombined (n = 20 + 11)")

ps_loop_net <- loop_net2+loop_net1+loop_net3+plot_layout(nrow=1) +
  plot_annotation(title="Discovery (N=35)", 
                  theme=theme(plot.title = element_text(face="bold",hjust=0.5)))
ps3 <- ps_loop_net + p7
ggsave("results/Train_Val_net.pdf", ps3, width=12, height=4)

# lasso ----
## epitranscriptome ----
model_nad_la <- glmNetFit(model_df1, model = "lasso",
                          y="age", seed=10, train.prop = 0.7, 
                          split = T, sample.group = meta.uni$age_group)

# get selected variables, intercept as the first element
sel_var_nad_la <- coef(model_nad_la$fit)[which(abs(as.vector(coef(model_nad_la$fit))) > 0),]

# LOO for prediction
nad_loo_la <- LOOP(model_df1, y="age", model = "lasso", 
                   alpha = model_nad_la$bestTune$alpha, 
                   lambda = model_nad_la$bestTune$lambda)
loop_la1 <- modelScatter(prediction = nad_loo_la$prediction, 
                         evalutaion = nad_loo_la$evaluation,
                         title = paste0("Epitranscriptome (n = ", length(sel_var_nad_la)-1,")"))

## transcriptome ----
model_expr_la <- glmNetFit(model_df2, model = "lasso",
                           y="age", seed=226, train.prop = 0.7, 
                           split = T, sample.group = meta.uni$age_group)

# get selected variables, intercept as the first element
sel_var_expr_la <- coef(model_expr_la$fit)[which(abs(as.vector(coef(model_expr_la$fit))) > 0),]

# LOO for prediction
expr_loo_la <- LOOP(model_df2, y="age", model = "lasso", 
                    alpha = model_expr_la$bestTune$alpha, 
                    lambda = model_expr_la$bestTune$lambda)
loop_la2 <- modelScatter(prediction = expr_loo_la$prediction, 
                         evalutaion = expr_loo_la$evaluation,
                         title = paste0("Transcriptome (n = ", length(sel_var_expr_la)-1,")"))

## epitranscriptome and transcriptome ----
expr_sel_la <- t(expr.norm[names(sel_var_expr_la)[-1],meta.uni$sample.id])
colnames(expr_sel_la) <- paste0(colnames(expr_sel_la), ".expr")
nad_sel_la <- t(ind.lfc.keep.fc2[names(sel_var_nad_la)[-1],meta.uni$sample.id])
colnames(nad_sel_la) <- paste0(colnames(nad_sel_la), ".nad")

# combined model 
model_df3_la <- data.frame(nad_sel_la, expr_sel_la, age=meta.uni$age)     
model_combined_la <- glmNetFit(model_df3_la, y="age", seed=11, 
                               model = "lasso",
                               train.prop = 0.7, split = T,
                               sample.group = meta.uni$age_group)
sel_var_combined_la <- coef(model_combined_la$fit)[which(abs(as.vector(coef(model_combined_la$fit))) > 0),]
nad_var_num_la <- length(grep("nad",names(sel_var_combined_la[-1])))
expr_var_num_la <- length(sel_var_combined_la)-nad_var_num_la-1

# LOO for prediction
combined_loo_la <- LOOP(model_df3_la, y="age", model = "lasso", 
                        alpha = model_combined_la$bestTune$alpha, 
                        lambda = model_combined_la$bestTune$lambda)
loop_la3 <- modelScatter(prediction = combined_loo_la$prediction, 
                         evalutaion = combined_loo_la$evaluation,
                         title = paste0("Combined (n = ",expr_var_num_la,"+",nad_var_num_la,")"))
## validation ----
valid.m.la <- validPred(nad.mat = ind.lfc.val, expr.mat = expr.val,
                        model = model_combined_la$fit, obs = metadata.uni.val$age)
p7.la <- modelScatter(valid.m.la$prediction, valid.m.la$evaluation, 
                      title="Validation (N=26)\nCombined (n = 9 + 8)")


ps_loop_la <- loop_la2+loop_la1+loop_la3+plot_layout(nrow=1) +
  plot_annotation(title="Discovery (N=35) -- Lasso", 
                  theme=theme(plot.title = element_text(face="bold",hjust=0.5)))

ps3.la <- ps_loop_la + p7.la
ggsave("results/reg/Train_Val_lasso.pdf", ps3.la, width=12, height=4)

# ridge ----
## epitranscriptome ----
model_nad_ri <- glmNetFit(model_df1, model = "ridge",
                          y="age", seed=100, train.prop = 0.7, 
                          split = T, sample.group = meta.uni$age_group)

# get selected variables, intercept as the first element
sel_var_nad_ri <- coef(model_nad_ri$fit)[which(abs(as.vector(coef(model_nad_ri$fit))) > 0),]

# LOO for prediction
nad_loo_ri <- LOOP(model_df1, y="age", model = "ridge", 
                   alpha = model_nad_ri$bestTune$alpha, 
                   lambda = model_nad_ri$bestTune$lambda)
loop_ri1 <- modelScatter(prediction = nad_loo_ri$prediction, 
                         evalutaion = nad_loo_ri$evaluation,
                         title = paste0("Epitranscriptome"))

## transcriptome ----
model_expr_ri <- glmNetFit(model_df2, model = "ridge",
                           y="age", seed=11, train.prop = 0.7, 
                           split = T, sample.group = meta.uni$age_group)

# get selected variables, intercept as the first element
sel_var_expr_ri <- coef(model_expr_ri$fit)[which(abs(as.vector(coef(model_expr_ri$fit))) > 0),]

# LOO for prediction
expr_loo_ri <- LOOP(model_df2, y="age", model = "ridge", 
                    alpha = model_expr_ri$bestTune$alpha, 
                    lambda = model_expr_ri$bestTune$lambda)
loop_ri2 <- modelScatter(prediction = expr_loo_ri$prediction, 
                         evalutaion = expr_loo_ri$evaluation,
                         title = paste0("Transcriptome"))

## epitranscriptome and transcriptome ----
expr_sel_ri <- t(expr.norm[names(sel_var_expr_ri)[-1], meta.uni$sample.id])
colnames(expr_sel_ri) <- paste0(colnames(expr_sel_ri), ".expr")
nad_sel_ri <- t(ind.lfc.keep.fc2[names(sel_var_nad_ri)[-1], meta.uni$sample.id])
colnames(nad_sel_ri) <- paste0(colnames(nad_sel_ri), ".nad")

# combined model 
model_df3_ri <- data.frame(nad_sel_ri, expr_sel_ri, age=meta.uni$age)
model_combined_ri <- glmNetFit(model_df3_ri, y="age", seed=10, 
                               model = "ridge",
                               train.prop = 0.7, split = T,
                               sample.group = meta.uni$age_group)

sel_var_combined_ri <- coef(model_combined_ri$fit)[which(abs(as.vector(coef(model_combined_ri$fit))) > 0),]
nad_var_num_ri <- length(grep("nad",names(sel_var_combined_ri[-1])))
expr_var_num_ri <- length(sel_var_combined_ri)-nad_var_num_ri-1

# LOO for prediction
combined_loo_ri <- LOOP(model_df3_ri, y="age", model = "ridge", 
                        alpha = model_combined_ri$bestTune$alpha, 
                        lambda = model_combined_ri$bestTune$lambda)
loop_ri3 <- modelScatter(prediction = combined_loo_ri$prediction, 
                         evalutaion = combined_loo_ri$evaluation,
                         title = paste0("Combined"))

## validation ----
valid.m.ri <- validPred(nad.mat = ind.lfc.val, expr.mat = expr.val,
                        model = model_combined_ri$fit, obs = metadata.uni.val$age)
p7.ri <- modelScatter(valid.m.ri$prediction, valid.m.ri$evaluation, 
                      title="Validation (N=26)\nCombined")

ps_loop_ri <- loop_ri2+loop_ri1+loop_ri3+plot_layout(nrow=1) +
  plot_annotation(title="Discovery (N=35) -- Ridge", 
                  theme=theme(plot.title = element_text(face="bold",hjust=0.5)))
ps3.ri <- ps_loop_ri + p7.ri
ggsave("results/reg/Train_Val_Ridge.pdf", ps3.ri, width=12, height=4)

# random forest regression model ----
## epitranscriptome ----
model_nad_rf <- RFRFit(model_df1,
                       y="age", seed=468, train.prop = 0.7, 
                       split = T, sample.group = meta.uni$age_group)
# LOO for prediction
nad_loo_rf <- LOOP(model_df1, y="age", model = "randomforest", 
                   mtry = model_nad_rf$bestTune$mtry)
loop_rf1 <- modelScatter(prediction = nad_loo_rf$prediction, 
                         evalutaion = nad_loo_rf$evaluation,
                         title = paste0("Epitranscriptome"))

## transcriptome ----
model_expr_rf <- RFRFit(model_df2,
                        y="age", seed=468, train.prop = 0.7, 
                        split = T, sample.group = meta.uni$age_group)
# LOO for prediction
expr_loo_rf <- LOOP(model_df2, y="age", model = "randomforest", 
                    mtry = model_expr_rf$bestTune$mtry)
loop_rf2 <- modelScatter(prediction = expr_loo_rf$prediction, 
                         evalutaion = expr_loo_rf$evaluation,
                         title = paste0("Transcriptome"))

## epitranscriptome and transcriptome ----
expr_sel_rf <- t(expr.norm[expr_cor, meta.uni$sample.id])
colnames(expr_sel_rf) <- paste0(colnames(expr_sel_rf), ".expr")
nad_sel_rf <- t(ind.lfc.keep.fc2[nad_cor, meta.uni$sample.id])
colnames(nad_sel_rf) <- paste0(colnames(nad_sel_rf), ".nad")

# combined model 
model_df3_rf <- data.frame(nad_sel_rf, expr_sel_rf, age=meta.uni$age)
model_combined_rf <- RFRFit(model_df3_rf,
                            y="age", seed=917, train.prop = 0.7, 
                            split = T, sample.group = meta.uni$age_group)
# LOO for prediction
combined_loo_rf <- LOOP(model_df3_rf, y="age", model = "randomforest", 
                        mtry = model_combined_rf$bestTune$mtry)
loop_rf3 <- modelScatter(prediction = combined_loo_rf$prediction, 
                         evalutaion = combined_loo_rf$evaluation,
                         title = paste0("Combined"))

## validation ----
valid.m.rf <- validPred(nad.mat = ind.lfc.val, expr.mat = expr.val,
                        model = model_combined_rf$fit, obs = metadata.uni.val$age,
                        model.coef = c(colnames(expr_sel_rf), colnames(nad_sel_rf)))

p7.rf <- modelScatter(valid.m.rf$prediction, valid.m.rf$evaluation, 
                      title="Validation (N=26)\nCombined")

ps_loop_rf <- loop_rf2+loop_rf1+loop_rf3+plot_layout(nrow=1) +
  plot_annotation(title="Discovery (N=35) -- Random Forest", 
                  theme=theme(plot.title = element_text(face="bold",hjust=0.5)))

ps3.rf <- ps_loop_rf + p7.rf
ggsave("results/reg/Train_Val_RandomForest_v4.pdf", ps3.rf, width=12, height=4)

# save model
save(model_nad, model_expr, model_combined,
     model_nad_la, model_expr_la, model_combined_la,
     model_nad_ri, model_expr_ri, model_combined_ri,
     model_nad_rf, model_expr_rf, model_combined_rf,
     file = "Data/AgeModel.RData")

# Coefficients used in age prediction model ----
## elastic net ----
### trans ----
net_var_df_tr <- data.frame(id = names(sel_var_expr),
                         coef = sel_var_expr) %>% 
  left_join(gene_anno, by=c("id"="ensembl_gene_id"))

### epitrans ----
net_var_df_epi <- data.frame(id = names(sel_var_nad),
                         coef = sel_var_nad) %>% 
  left_join(gene_anno, by=c("id"="ensembl_gene_id"))

### combined ----
net_var_df <- data.frame(id = names(sel_var_combined[-1]),
                         coef = sel_var_combined[-1])
net_var_df$type <- gsub(".*\\.", "", net_var_df$id)
net_var_df$id <- gsub("\\..*", "", net_var_df$id)
net_var_df <- net_var_df %>% 
  left_join(gene_anno, by=c("id"="ensembl_gene_id"))

# ENSG00000250461: FAM133BPS
net_var_df$external_gene_name[net_var_df$id == "ENSG00000250461"] <- "FAM133BP"

# ENSG00000260257: NOVELGENE
net_var_df$external_gene_name[net_var_df$id == "ENSG00000260257"] <- "NOVELGENE"

# ENSG00000223612: ECOPP
net_var_df$external_gene_name[net_var_df$id == "ENSG00000223612"] <- "ECOPP"

# for GPATCH4 NAD-RNA: GPATCH4.nad
net_var_df$external_gene_name[net_var_df$external_gene_name == "GPATCH4" & net_var_df$type == "nad"] <- "GPATCH4.nad"

net_var_df$external_gene_name <- factor(net_var_df$external_gene_name, 
                                        levels=net_var_df$external_gene_name[order(net_var_df$coef)])

cp1 <- net_var_df %>% 
  ggplot(aes(external_gene_name, y=1)) +
  geom_point(aes(shape=type, color=coef), size=10) +
  coord_polar() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color="black", vjust=1),
        legend.position = "bottom") +
  ylim(c(0, 1)) +
  scale_color_gradientn(colors=rev(paint_palette("Pearlgirl", 100, "continuous"))) +
  labs(x='', y='', shape='', color="Coefficients") 
ggsave("results/EL_coef.pdf", cp1, width = 8, height = 8)

# draw top 3 pos. and neg. coef. genes
it.genes2 <- rbind(slice_max(net_var_df, order_by = coef, n=3), 
                   slice_min(net_var_df, order_by = coef, n=3))
pls1 <- list()
for (i in 1:nrow(it.genes2)) {
  if (it.genes2[i,]$type == "expr") {
    df1 <- expr.norm
  } else {
    df1 <- ind.lfc.keep.fc2
  }
  
  p1 <- data.frame(value=df1[it.genes2[i,]$id,],
                   id=colnames(df1)) %>% 
    mutate(group=factor(gsub("\\d+","",id), levels=c("Y", "M", "O"))) %>% 
    ggplot(aes(group, value)) +
    geom_boxplot(aes(fill=group)) + 
    ggsignif::geom_signif(comparisons = list(c("Y","M"),c("Y","O"),c("M","O")),
                          tip_length = 0, step_increase = 0.1, parse = T, 
                          map_signif_level = function(p) sprintf("italic(P) == %.2g", p)) +
    theme_classic() +
    theme(axis.text.y = element_text(color='black'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(face="bold", hjust=0.5),
          legend.position = "bottom") +
    scale_fill_manual(values=c("#edd0ca","#aa688f","#2d1d3d"),
                      labels=c("Young","Mid","Old")) +
    labs(x="", fill="", subtitle = it.genes2[i,]$type, 
         title = it.genes2[i,]$external_gene_name)
  pls1[[i]] <- p1
}

ps1 <- wrap_plots(pls1) + 
  plot_layout(nrow=2, guides = "collect") &
  theme(legend.position='bottom')

cp1 + ps1 + plot_layout(widths = c(1.5, 2))
ggsave("results/reg/EL_coef_top6.pdf", width = 16, height = 8)