# Validation of aging clock
library(enONE)
library(tidyverse)
library(patchwork)
library(caret)
library(glmnet)

# Custom Function #
# glm regression with elastic net
glmNetFit <- function(data, x=NULL, y) {
  pred.var <- model.matrix(formula(paste(y, "~ .")), data=data)[,-1]
  resp.var <- data[, y]
  
  model_train <- caret::train(formula(paste(y, "~ .")),
                              data = data, method = "glmnet",
                              trControl = caret::trainControl(method="LOOCV")
  )
  # fit final model
  net.model <- glmnet(pred.var, resp.var, 
                      family="gaussian",
                      alpha = model_train$bestTune$alpha, 
                      lambda = model_train$bestTune$lambda)
  
  p.net <- eval.mat <- NULL
  # get prediction 
  p.net <- subset(model_train$pred, alpha == model_train$bestTune$alpha & lambda == model_train$bestTune$lambda)
  
  # evaluation
  # model performance evaluation
  eval.mat <- evalMetrics(obs = p.net$obs, pred = p.net$pred)
  
  return(list(
    fit = net.model,
    prediction = p.net,
    evaluation = eval.mat,
    bestTune = model_train$bestTune
  ))
}

# calculate evaluation metrics from observation and prediction
evalMetrics <- function(obs, pred) {
  # delta age 
  dage <- pred - obs
  
  # evaluation
  # model performance evaluation
  eval.mat <- data.frame(
    RMSE = RMSE(pred, obs),
    Rsquare = R2(pred, obs),
    MAD = median(abs(dage)),
    Cor = cor.test(pred, obs)$estimate,
    Cor.p = cor.test(pred, obs)$p.value
  )
}

# draw scatter plot from `glmNetFit` output summarizing model performance
modelScatter <- function(model.list, title=NULL) {
  
  p.format <-  rstatix::p_format(model.list$evaluation$Cor.p, 
                                 digits = 3, leading.zero = FALSE,
                                 trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16)
  
  # associated model evaluation metrics 
  lab <- paste(paste0("MAE=", round(model.list$evaluation$MAE,3)), 
               paste0("cor=", round(model.list$evaluation$Cor,3)),
               p.format, sep="; ")
  
  p1 <- ggplot(model.list$prediction, aes(obs, pred)) +
    geom_segment(x=min(model.list$prediction$obs), y=min(model.list$prediction$obs), 
                 xend=max(model.list$prediction$obs), yend=max(model.list$prediction$obs), 
                 lty="dashed") +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic() + 
    theme(axis.text = element_text(color="black")) +
    labs(x="Chronological Age", y="Predicted Age", 
         subtitle = lab, title=title)
  
  return(p1)
}
# predict with validation set
validPred <- function(nad.mat=NULL, expr.mat=NULL, model, obs) {
  coef_model <- rownames(coef(model))[-1]
  
  if (!is.null(nad.mat) & !is.null(expr.mat)) {
    nad.sel <- t(nad.mat[gsub("\\..*","",grep("nad", coef_model, value=T)),])
    colnames(nad.sel) <- paste0(colnames(nad.sel), ".nad")
    
    expr.sel <- t(expr.mat[gsub("\\..*","",grep("expr", coef_model, value=T)),])
    colnames(expr.sel) <- paste0(colnames(expr.sel), ".expr")
    
    test.df <- data.frame(nad.sel, expr.sel, age=obs)
    
  } else if (is.null(expr.mat) & !is.null(nad.mat)) {
    nad.sel <- t(nad.mat[coef_model,])
    test.df <- data.frame(nad.sel, age=obs)
    
  } else if (!is.null(expr.mat) & is.null(nad.mat)) {
    expr.sel <- t(expr.mat[coef_model,])
    test.df <- data.frame(expr.sel, age=obs)
    
  }
  
  x.test <- model.matrix(age ~., test.df)[,-1]
  
  val.p <- as.vector(predict(model, x.test))
  
  pred.df <- data.frame(obs = obs, pred = val.p)
  # calculate evaluation metrics
  eval.mat <- evalMetrics(obs=pred.df$obs, pred=pred.df$pred)
  return(list(
    fit = model,
    prediction = pred.df,
    evaluation = eval.mat
  ))
  
}

# load data, including Enone, top.norm.data, and res.sig.ls (NAD-RNAs)
load("DATA/Enone.RData")

# Enone workflow for validation cohort
## set up ----
## load data, set spike-in prefix, set sample control index
# In Data folder 
metadata <- read.csv("Data/metadata.csv")
counts_df <- read.csv("Data/Counts.csv", row.names = 1)

# use validation cohort
metadata.val <- metadata[metadata$cohort.group == "Validation",]
counts.val <- counts_df[, metadata.val$id]

metadata.uni.val <- metadata.val[!duplicated(metadata.val$sample.id), ]
# prefix of Drosophila spike-in gene id
spikeInPrefix <- "^FB"

## filtering ----
counts.val.keep <- counts.val[rownames(Enone),]

## create Enone ----
Enone.val <- createEnone(counts.val.keep, 
                     bio.group = metadata.val$condition,
                     enrich.group = metadata.val$Assay,
                     batch.group = metadata.val$batch,
                     spike.in.prefix = spikeInPrefix,
                     synthetic.id = c("Syn1","Syn2"),
                     input.id = "Input",
                     enrich.id = "Enrich")
## run ---- 
Enone.val <- enONE(Enone.val, 
               ruv.norm = TRUE, ruv.k = 5, 
               eval.pc.n = 5, eval.pam.k = 2:12,
               return.norm = TRUE)

## check performance
enScore <- getScore(Enone.val)

# get normalized counts
counts.norm.val <- Counts(Enone.val, "sample", "DESeq_RUVs_k5")
# individual enrichment profiles ----
counts.norm.log.val <- log2(counts.norm.val + 1)
enrich_idx <- CreateGroupMatrix(Enone.val$enrich)

# log
ind.lfc.val <- counts.norm.log.val[, enrich_idx[1,]] - counts.norm.log.val[, enrich_idx[2,]]
colnames(ind.lfc.val) <- metadata.uni.val$sample.id

# expression
expr.val <- counts.norm.log.val[,enrich_idx[2,]]
colnames(expr.val) <- metadata.uni.val$sample.id

# validation ----
# load model, including: 
# model_expr, model_expr, model_combined_final, 
# sel_var_nad, sel_var_expr, sel_var_combined_final
load("DATA/AgeModel.RData")

# combined model: 34 = 20 + 14
valid.m <- validPred(nad.mat = ind.lfc.val, expr.mat = expr.val, 
                      model = model_combined_final$fit, obs = metadata.uni.val$age)
p1 <- modelScatter(valid.m, title="Validation (N=8)\nCombined model (n = 20+14)")
ggsave("results/Val_combined_model.pdf", p1, width=6, height=5)
