# Identification of NAD-RNA biomarkers
# Construct a logistic regression model with an elastic net penalty based on the scaled NAD modification levels
library(enONE)
library(tidyverse)
library(paintingr)
library(caret)
library(glmnet)
library(pROC)

# load NAD-RNA data from results of 03_enONE.R
load('Data/Enone.RData')
# load related data from results of 05-1_NAD-RNA_dynamics.R
load('Data/StatDat.RData')

# computing error metrics given a confusion matrix
err_metric <- function(CM) {
  TN = CM[1,1]
  TP = CM[2,2]
  FP = CM[1,2]
  FN = CM[2,1]
  sensitivity = (TP)/(TP+FN)
  specificity = (TN)/(TN+FP)
  accuracy_model  = (TP+TN)/(TP+TN+FP+FN)
  print(paste("Sensitivity value of the model: ",round(sensitivity,2)))
  print(paste("Specificity value of the model: ",round(specificity,2)))
  print(paste("Accuracy of the model: ",round(accuracy_model,2)))
  # return error metrics
  err_met <- data.frame(metrics = c("TP","TN","FP","FN",,
                                    "Sensitivity","Specificity","Accuracy"),
                        values = c(TP,TN,FP,FN,
                                   sensitivity,specificity,accuracy_model)
  )
  return(err_met)
}

# differential modified NAD-RNA ----
# s3: up-regulated NAD-RNAs upon NMN; s4: down-regulated NAD-RNAs upon NMN
# get differential modified NAD-RNAs table of modification levels
nad_df_sub <- subset(nad_df, GeneID %in% c(s3,s4))
nad_df_sub <- nad_df_sub %>% dplyr::group_by(Group) %>% distinct(GeneID, .keep_all = TRUE) %>% as.data.frame()

bvp1 <- BetweenStatPlot(nad_df_sub, x='Group', y='logFC', color='Group',
                comparisons = list(c('D0_ctrl','D14_ctrl'),c('D0_nmn','D14_nmn')),
                palette = c('#56565C','#A9A89F','#176AA6','#4F99B4'),
                step.increase = 0.1) + labs(y='Log2 Fold Change')

# get differential modified NAD-RNAs table of experssion levels
expr_df_sub <- subset(expr_df, GeneID %in% c(s3,s4))
expr_df_sub <- expr_df_sub %>% dplyr::group_by(Group) %>% distinct(GeneID, .keep_all = TRUE) %>% as.data.frame()

bvp2 <- BetweenStatPlot(expr_df_sub, x='Group', y='logexpr', color='Group',
                comparisons = list(c('D0.ctrl','D14.ctrl'),c('D0.nmn','D14.nmn')),
                palette = c('#56565C','#A9A89F','#176AA6','#4F99B4'),
                step.increase = 0.3) + labs(y='Log2 Normalized counts')

# fig. 5A
# plot expression and enrichment of differential modified NAD-RNA
bvps1 <- bvp2+bvp1
ggsave('DynamicNAD_FC_Expr.pdf', bvps1, width=6, height=4)

# get significant enrichment
res.sig.ls <- getEnrichment(Enone, slot='sample', filter=T)
res.sig.ls <- lapply(res.sig.ls, function(x) subset(x, logCPM>1))

# logistic regression with elastic net penalty ----
# 'a' encoded without NMN exposure; 'b' encoded with NMN exposure
nmn.response <- c(rep('a',13), rep('b',5))

# split data into training and test data, take 70% of the data as training set 
sample_size <- floor(0.7*ncol(ind.fc)) # ind.fc: sample-wise log fold-change matrix

# set seed ensure the sample function produce same random numbers
set.seed(1212)
training_index <- sample(seq_len(ncol(ind.fc)), size = sample_size)

## pre-selecting features ----
# use those differential modified genes with logCPM >= 2.5
cand.set <- c(s3, s4) # size of (s3+s4) = 392
cand.set <- cand.set[cand.set %in% subset(res.sig.ls$D0_ctrl.Enrich_D0_ctrl.Input, logCPM>=2.5)$GeneID]
# size of cand.set = 151

## model training ----
# prepare data
model_df <- data.frame(t(ind.fc[cand.set,]), resp=nmn.response)
model_df$resp <- as.factor(model_df$resp)

# training data (70%)
mkr.train <- model_df[training_index, ]
x <- model.matrix(resp~., data=mkr.train)[,-1]
y <- as.numeric(mkr.train$resp)-1

# test data (30%)
mkr.test <- model_df[-training_index, ]
x.test <- model.matrix(resp~., data=mkr.test)[,-1]

# performing stratified CV due to insufficient samples 
# create sample id for each fold
set.seed(1104)
fold.id <- createFolds(y, k=5, returnTrain = TRUE)

# determine the hyperparameters: alpha and lambda 
model.net <- caret::train(
  resp ~., data = mkr.train, method='glmnet', metric ='ROC',
  trControl = trainControl("cv", number=5, index=fold.id, classProbs = TRUE)
)
model.net$bestTune

# build final model   
log.net.model <- glmnet(x, y, alpha = model.net$bestTune$alpha, family = 'binomial',
                        lambda = model.net$bestTune$lambda)
table(abs(as.vector(coef(log.net.model))) > 0)

# get selected variables, excluding intercept (first element)
sel.var <- coef(log.net.model)[which(abs(as.vector(coef(log.net.model))) > 0),][-1]

## prediction ----
p.ln <- predict(log.net.model, x.test)
predicted.classes <- ifelse(p.ln > 0.5, 'b', 'a') # probability threshold = 0.5

# model performance evaluation ----
# get confusion matrix
observed.classes <- mkr.test$resp
conf.mat <- table(observed.classes,predicted.classes)

# comput error metrics
err.met <- err_metric(conf.mat)

# 90% CI of accuracy
# success: TP+TN = 5; cases TP+TN+FP+FN = 6
ratesci::jeffreysci(5, 6, level = 0.9)

# ROC 
roc_score <- pROC::roc(mkr.test$resp, predicted.classes[,1])
# 90% bootstrap CI of AUC
tmp <- pROC::ci.auc(mkr.test$resp, predicted.classes[,1], 
                   conf.level=0.90, method='bootstrap', boot.n=2000)

# ROC curve
roc.df <- data.frame(x=1-rev(roc_score$specificities), y=rev(roc_score$sensitivities))
rocp1 <- ggplot(roc.df, aes(x, y)) +
  geom_line(aes(group=1), size=1, color='#19669d') +
  geom_abline(slope = 1, lty = 'dashed', color='grey40') +
  annotate('text', x=0.8, y=0.1, label='AUC=0.75 [0.5-1]') +
  annotate('text', x=0.7, y=0.2, label='Accuaracy=0.83 [0.50-0.96]') +
  theme_classic() +
  theme(axis.text = element_text(color='black')) +
  labs(x='1-Specificity', y='Sensitivity')
ggsave('ROC_log_net.pdf', rocp1, width=4, height=3)

# save model
save(ind.fc, log.net.model, file='data/biomrk.RData')

# hierarchical clustering based on 4 NAD-RNA biomarkers ----
# get symbol of biomarkers
sel.symbol <- setNames(subset(gene_id, ensembl_gene_id %in% names(sel.var))$external_gene_name,
                       nm=subset(gene_id, ensembl_gene_id %in% names(sel.var))$ensembl_gene_id)
sel.symbol <- sel.symbol[!duplicated(names(sel.symbol))]
sel.symbol <- sel.symbol[names(sel.var)]
names(sel.symbol) <- names(sel.var)

# missing annotation
sel.symbol['ENSCAFG00845015316'] <- 'NovelGene'

# scaled NAD-RNA modification levels of biomarkers
sel.fc <- t(scale(t(ind.fc[names(sel.var),])))
rownames(sel.fc) <- sel.symbol

# palette for heatmap
col_fun <- circlize::colorRamp2(breaks=c(-2,-1,0,1,2), 
                                colors=c("#4575B4","#ABD9E9","#FFFFFF","#FDAE61","#dc443c"))
# column annotation
col_anno_all <- HeatmapAnnotation(group = c(rep(c('D0_Water', 'D14_Water'), each=4), 
                                           rep(c('D0_NMN', 'D14_NMN'), each=5)),
                                  col = list(group = c('D0_Water'='#D3D4D0','D14_Water'='#040000',
                                                       'D0_NMN'='#9595C5','D14_NMN'='#4B4C6B')))
# heatmap of scaled NAD-RNA modification levels 
hm1 <- Heatmap(sel.fc, name = 'z-score', col = col_fun, top_annotation = col_anno_all)

# scaled expression levels of biomarkers
ind.expr <- counts_norm[,enrich_idx[2,]]
colnames(ind.expr) <- unique(gsub('(Input.)|(Enrich.)','',samples_name))
sel.expr <- t(scale(t(ind.expr[names(sel.var),])))
rownames(sel.expr) <- sel.symbol
hm2 <- Heatmap(sel.expr, name = 'z-score', col = col_fun, top_annotation = col_anno_all)
# draw heatmaps
pdf('NAD_biomrks_hm.pdf', width=4, height=3)
draw(hm1)
dev.off()

pdf('Expr_biomrks_hm.pdf', width=4, height=3)
draw(hm2)
dev.off()

# create matrix for mis-classified samples of hierarchical clustering 
# In each vector, first element: FP and second element: FN
# In the list, first element: for enrichment hc and second element: for expression hc
n_mis_hc_ls <- list(c(0,1), c(1,3))
n_samples_case <- c(13, 5)
# confusion matrix
conf_mat_ls <- lapply(n_mis_hc_ls, function(x) {
  rbind(
    c(n_samples_case[1] - x[1], x[1]),
    c(x[2], n_samples_case[2] - x[2])
  )
})

hc_stats_ls <- lapply(conf_mat_ls, err_metric)

# performance evaluation for hclustering
# For epitranscriptome
# 90% CI of sensitivity
# success: TP = 4; cases: TP+FN = 5
ratesci::jeffreysci(4, 5, level = 0.9)

# 90% CI of specificity
# success: TN = 13; cases TN+FP = 13
ratesci::jeffreysci(13, 13, level = 0.9)

# For transcriptome
# 90% CI of sensitivity
# success: TP = 2; cases: TP+FN =5
ratesci::jeffreysci(2, 5, level = 0.9)

# 90% CI of specificity
# success: TN = 12; cases TN+FP = 13
ratesci::jeffreysci(12, 13, level = 0.9)

# logistic regression on each biomarker ----
# training index (70%)
set.seed(1106)
training_index_sg <- sample(seq_len(ncol(ind.fc)), size = sample_size)

# prepare data
# 0 encoded  without NMN exposure; 1 encoded with NMN exposure.
resp <- c(rep(0,13), rep(1,5))
resp.train <- as.factor(resp[training_index_sg])
resp.test <- as.factor(resp[-training_index_sg])

# initialize data 
model.ls <- list() # for stroing models
roc_df_sel <- data.frame() # for ROC metrics
# model training and prediction for each biomarker
for (i in names(sel.var)) {
  # model training 
  x.train <- as.vector(scale(ind.fc[i, training_index_sg]))
  x.test <- as.vector(scale(ind.fc[i, -training_index_sg]))
  log.lm1 <- glm(resp.train ~ x.train, family = 'binomial')
  model.ls[[ sel.symbol[i] ]] <- log.lm1
  p.model <- summary(log.lm1)[["coefficients"]][2,"Pr(>|z|)"]

  # prediction
  ## with test set
  p.n <- predict(log.lm1, data.frame('x.train'=x.test), type='response')
  p.class <- ifelse(p.n > 0.5, 1, 0)
  
  # ROC
  roc_score_i <- pROC::roc(resp.test, p.class)
  # 90% bootstrap CI of AUC
  roc_ci <- pROC::ci.auc(resp.test, p.class, conf.level=0.90, 
                        method='bootstrap', boot.n=2000)
  roc_ci <- round(roc_ci,2)
  roc_label <- paste0(sel.symbol[i]," (AUC = ", roc_ci[2], " [", roc_ci[1],"-", roc_ci[3], "])")
  roc_df_i <- data.frame(x=1-rev(roc_score_i$specificities), 
                         y=rev(roc_score_i$sensitivities),
                         symbol=sel.symbol[i],
                         label=roc_label) 

  # each row is a biomarker and its model performance
  roc_df_sel <- rbind(roc_df_sel, roc_df_i) 

}

# roc curve
roc_df_sel$symbol <- factor(roc_df_sel$symbol, levels = unique(roc_df_sel$symbol))
rocps1 <- ggplot(roc_df_sel, aes(x, y, group=symbol, color=symbol)) +
  geom_line(size=1) +
  geom_abline(slope = 1, lty = 'dashed', color='grey40') +
  theme_classic() +
  theme(axis.text = element_text(color='black')) +
  scale_color_manual(values=paint_palette('Twilight',4,'continuous'),
                     label=unique(roc_df_sel$label)) +
  labs(x='1-Specificity', y='Sensitivity', color='')

# transform to long format
sel.fc.long <- as.data.frame(sel.fc) %>% 
  rownames_to_column('symbol') %>% 
  pivot_longer(cols = -symbol,
               names_to = 'id',
               values_to = 'z') %>% 
  dplyr::mutate(Group=factor(gsub('\\..*', '', id), levels=c('D0_ctrl','D14_ctrl','D0_nmn','D14_nmn')),
                Exp=if_else(Group %in% 'D14_nmn', 'Exp', 'NoExp'))

sel.fc.long$Exp <- factor(sel.fc.long$Exp, levels=c('NoExp', 'Exp'))

# boxplot of selected NAD-RNA modificaion levels in z-score
bp1 <- ggplot(sel.fc.long, aes(Exp, z)) +
  geom_boxplot(aes(color=Exp)) +
  ggsignif::geom_signif(comparisons = list(c('Exp','NoExp')), 
                        test = 'wilcox.test', tip_length = 0) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color='black'),
        strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, "lines")) +
  facet_wrap(~symbol, nrow=1, strip.position = 'bottom', ) +
  scale_color_manual(values = c('black', '#e7aa0d')) +
  labs(x='', y='z-score', color='')

bp1 + rocps1 
ggsave('Biomrks_reg_enrich.pdf', width = 10, height = 4)



