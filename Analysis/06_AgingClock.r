# Age prediction model 
# based on discovery cohort
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

# load data, including nad_id, ind.lfc.keep.fc2, ind.scale.fc2, expr.norm, expr.scale
load("DATA/NADRNA_profiles.RData")

# In Data folder 
metadata <- read.csv("Data/metadata.csv")

# only use discovery cohort
metadata <- metadata[metadata$cohort.group == "Discovery",]
metadata.uni <- metadata[!duplicated(metadata$sample.id),]

# modeling ----
## model based on epitranscriptome ----
# prepare data
model_df1 <- data.frame(t(ind.lfc.keep.fc2), age=metadata.uni$age)

# fit model
model_nad <- glmNetFit(model_df1, y="age")

# get selected variables, intercept as the first element
sel_var_nad <- coef(model_nad$fit)[which(abs(as.vector(coef(model_nad$fit))) > 0),]

## model based on transcriptome ----
# prepare data
# compute cv2 to define highly variable genes
expr_avg <- rowMeans(2^expr.norm)
expr_var <- MatrixGenerics::rowVars(2^expr.norm)
expr_cv2 <- 2^(log2(expr_var) - log2(expr_avg^2))
expr_cv2 <- sort(expr_cv2, decreasing = TRUE)

# use top 2000 variable genes
model_df2 <- data.frame(t(head(expr.norm[names(expr_cv2),], n=2000)), age=metadata.uni$age)

# fit model
model_expr <- glmNetFit(model_df2, y="age")

# get selected variables, intercept as the first element
sel_var_expr <- coef(model_expr$fit)[which(abs(as.vector(coef(model.expr$fit))) > 0),]

# draw scatter plot
p1 <- modelScatter(model.expr, title="Transcriptome (n = 45)")
p2 <- modelScatter(model.nad, title="NAD-modified Epitranscriptome (n = 44)")

## model based on epitranscriptome and transcriptome ----
expr_sel <- t(expr.norm[names(sel_var_expr)[-1],])
colnames(expr_sel) <- paste0(colnames(expr_sel), ".expr")
nad_sel <- t(ind.lfc.keep.fc2[names(sel_var_nad)[-1],])
colnames(nad_sel) <- paste0(colnames(nad_sel), ".nad")

# first combined model 
model_df3 <- data.frame(nad_sel, expr_sel, age=metadata.uni$age)
model_combined <- glmNetFit(model_df3, y="age")
sel_var_combined <- coef(model_combined$fit)[which(abs(as.vector(coef(model_combined$fit))) > 0),]

# post selection of features to reduce features
# randomly draw features from combined model
## use seed 280 and size 40
set.seed(280) 
sample_vars <- sample(names(sel_var_combined)[-1], size=40)
sample_df <- model_df3[,c(sample_vars,"age")]
model_combined_final <- glmNetFit(sample_df, y="age")
sel_var_combined_final <- coef(model_combined_final$fit)[which(abs(as.vector(coef(model_combined_final$fit))) > 0),]

p3 <- modelScatter(model_combined_final, title="Combined (n = 20+14)")

ps1 <- p1 + p2 + p3 + plot_layout(nrow=1) +
  plot_annotation(title="Discovery (N=61)", 
                  theme=theme(plot.title = element_text(face="bold",hjust=0.5)))
ggsave("results/AgingClockDiscovery.pdf", ps1, width=14, height=5)

# save
save(model_expr, model_expr, model_combined_final, 
     sel_var_nad, sel_var_expr, sel_var_combined_final,
     file="DATA/AgeModel.RData")