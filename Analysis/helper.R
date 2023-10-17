# Custom Function #
# glm regression 
glmNetFit <- function(data, x=NULL, y, 
                      split=TRUE,
                      model = c("glmnet", "lasso", "ridge"),
                      seed=123, train.prop=0.7, 
                      alpha=NULL, lambda=NULL,
                      sample.group=NULL) {
  
  if (split) {
    # split dataset to training and test set
    # prepare data
    set.seed(seed)
    if (!is.null(sample.group)) {
      tr_idx <- createDataPartition(sample.group, p=train.prop, list=F)  
    } else (
      tr_idx <- createDataPartition(1:ncol(data), p=train.prop, list=F)  
    )
    
    tr_set <- data[tr_idx,] # training set
    test_set <- data[-tr_idx,] # test set
    
    pred.var <- model.matrix(formula(paste(y, "~ .")), data=tr_set)[,-1]
    resp.var <- tr_set[, y]
    
    # which model to fit
    if (model == "glmnet") {
      # Define parameter grid to search over
      tuneGrid <- expand.grid(alpha = c(0.01, 0.05, 0.1, 0.5, 1),
                              lambda = seq(0.1, 10, length.out = 10))
      
      
    } else if (model == "lasso") {
      # Apply Lasso regression with LOOCV to the training set
      # Define parameter grid to search over
      tuneGrid <- expand.grid(alpha = 1,
                              lambda = seq(0.1, 10, length.out = 10))
      
      
    } else if (model == "ridge") {
      # Apply ridge regression with LOOCV to the training set
      # Define parameter grid to search over
      tuneGrid <- expand.grid(alpha = 0,
                              lambda = seq(0.1, 10, length.out = 10))
    }
    
    # Train model 
    model_train <- caret::train(formula(paste(y, "~ .")),
                                data = tr_set, method = "glmnet",
                                tuneGrid = tuneGrid,
                                trControl = caret::trainControl(method="LOOCV")
                                
    )
    
    if (is.null(alpha) & is.null(lambda)) {
      alpha <- model_train$bestTune$alpha
      lambda <- model_train$bestTune$lambda
    }
    
    
    # Build final model
    net.model <- glmnet(pred.var, resp.var, 
                        family="gaussian",
                        alpha = alpha, 
                        lambda = lambda)
    
    # Make predictions on testing data
    test_predictions <- predict(net.model, as.matrix(test_set[,-ncol(test_set)]))
    
    p.net <- eval.mat <- NULL
    
    # get prediction 
    # for training set
    p.net.train <- subset(model_train$pred, alpha == model_train$bestTune$alpha & lambda == model_train$bestTune$lambda)
    colnames(p.net.train) <- c("obs", "pred")
    
    # for test set
    p.net.test <- data.frame(obs = test_set$age, pred = test_predictions)
    colnames(p.net.test) <- c("obs", "pred")
    
    p.net <- list("train" = p.net.train, "test" = p.net.test)
    
    # evaluation
    # model performance evaluation
    eval.mat.train <- evalMetrics(obs = p.net.train$obs, pred = p.net.train$pred)
    eval.mat.test <- evalMetrics(obs = p.net.test$obs, pred = p.net.test$pred)
    
    eval.mat <- list("train" = eval.mat.train, "test" = eval.mat.test)
    
  } else {
    # use all data for training
    pred.var <- model.matrix(formula(paste(y, "~ .")), data=data)[,-1]
    resp.var <- data[, y]
    
    # which model to fit
    if (model == "glmnet") {
      # Define parameter grid to search over
      tuneGrid <- expand.grid(alpha = c(0.01, 0.05, 0.1, 0.5, 1),
                              lambda = seq(0.1, 10, length.out = 10))
      
      
    } else if (model == "lasso") {
      # Apply Lasso regression with LOOCV to the training set
      # Define parameter grid to search over
      tuneGrid <- expand.grid(alpha = 1,
                              lambda = seq(0.1, 10, length.out = 10))
      
      
    } else if (model == "ridge") {
      # Apply ridge regression with LOOCV to the training set
      # Define parameter grid to search over
      tuneGrid <- expand.grid(alpha = 0,
                              lambda = seq(0.1, 10, length.out = 10))
    }
    
    # Train model 
    model_train <- caret::train(formula(paste(y, "~ .")),
                                data = data, method = "glmnet",
                                tuneGrid = tuneGrid,
                                trControl = caret::trainControl(method="LOOCV")
                                
    )
    
    if (is.null(alpha) & is.null(lambda)) {
      alpha <- model_train$bestTune$alpha
      lambda <- model_train$bestTune$lambda
    }
    
    # Build final model
    net.model <- glmnet(pred.var, resp.var, 
                        family="gaussian",
                        alpha = alpha, 
                        lambda = lambda)
    
    p.net <- eval.mat <- NULL
    # get prediction
    p.net <- list("train"=subset(model_train$pred, alpha == model_train$bestTune$alpha & lambda == model_train$bestTune$lambda))
    
    # evaluation
    # model performance evaluation
    eval.mat <- list("train"=evalMetrics(obs = p.net[["train"]]$obs, pred = p.net[["train"]]$pred))
  }
  
  
  # save best tune hyperparameters
  bt <- data.frame(alpha=alpha, lambda=lambda)
  
  return(list(
    fit = net.model,
    prediction = p.net,
    evaluation = eval.mat,
    bestTune = bt
  ))
  
}

# Fit random forest regression (RFR) model
RFRFit <- function(data, x=NULL, y, 
                   split=TRUE, seed=123, train.prop=0.7, 
                   sample.group=NULL) {
  if (split) {
    # split dataset to training and test set
    # prepare data
    set.seed(seed)
    if (!is.null(sample.group)) {
      tr_idx <- createDataPartition(sample.group, p=train.prop, list=F)  
    } else (
      tr_idx <- createDataPartition(1:ncol(data), p=train.prop, list=F)  
    )
    
    tr_set <- data[tr_idx,] # training set
    test_set <- data[-tr_idx,] # test set
    
    pred.var <- model.matrix(formula(paste(y, "~ .")), data=tr_set)[,-1]
    resp.var <- tr_set[, y]
    
  } else {
    # use all data for training
    pred.var <- model.matrix(formula(paste(y, "~ .")), data=data)[,-1]
    resp.var <- data[, y]
  }
  
  # Define a sequence of mtry values to search over
  mtry.seq <- floor(seq(10, 100, length.out = 15))
  
  # Create an empty vector to store the LOOCV errors
  cv.errors <- rep(0,length(mtry.seq))
  
  # Create an empty list to store the LOOCV predictions
  cv.p.ls <- list()
  
  set.seed(seed)
  # Perform LOOCV-based hyperparameter tuning
  idx <- 1:length(resp.var)
  for (i in 1:length(mtry.seq)) {
    # Create an empty vector to store the LOOCV predictions
    cv.pred <- rep(0,length(resp.var))
    
    # Perform LOOCV with the current mtry value
    for (j in idx) {
      # Split the data into training and test sets
      train.index <- idx[-j]
      test.index <- idx[j]
      train.data <- pred.var[train.index,]
      train.labels <- resp.var[train.index]
      test.data <- pred.var[test.index,]
      
      # Fit the random forest model using the training set
      rf.model <- randomForest(train.data, train.labels, ntree=1000, mtry=mtry.seq[i])
      
      # Make predictions using the test set
      cv.pred[j] <- predict(rf.model, newdata=test.data)
    }
    cv.p.ls[[i]] <- cv.pred
    # Calculate the LOOCV error for the current mtry value
    cv.errors[i] <- mean((resp.var - cv.pred)^2)
  }
  
  # Get the optimal mtry value
  opt.mtry <- mtry.seq[which.min(cv.errors)]
  
  # Fit the final random forest model using the entire discovery cohort with the optimal mtry value
  rf.final <- randomForest(pred.var, resp.var, ntree=1000, mtry=opt.mtry)
  
  # get prediction
  p.net.train <- p.net.test <- eval.train <- eval.test <- NULL
  
  # for training set
  p.net.train <- data.frame(obs = resp.var, pred = cv.p.ls[[which.min(cv.errors)]])
  colnames(p.net.train) <- c("obs", "pred")
  
  # evaluation metrics
  eval.train <- evalMetrics(p.net.train$obs, p.net.train$pred) 
  
  if (split) {
    # for test set
    test_predictions <- predict(rf.final, test_set[, -ncol(test_set)])
    p.net.test <- data.frame(obs = test_set$age, pred = test_predictions)
    colnames(p.net.test) <- c("obs", "pred")
    
    # evaluation metrics
    eval.test <- evalMetrics(p.net.test$obs, p.net.test$pred) 
  }
  
  p.net <- list("train" = p.net.train, "test" = p.net.test)
  eval.mat <- list("train" = eval.train, "test" = eval.test)
  bt <- data.frame(mtry=opt.mtry)
  
  return(list(
    fit = rf.final,
    prediction = p.net,
    evaluation = eval.mat,
    bestTune = bt
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
    MAE = median(abs(dage)),
    Cor = cor.test(pred, obs)$estimate,
    Cor.p = cor.test(pred, obs)$p.value
  )
}

# draw scatter plot from `glmNetFit` output summarizing model performance
modelScatter <- function(prediction, evalutaion, title=NULL) {
  
  p.format <- rstatix::p_format(evalutaion$Cor.p, 
                                digits = 3, leading.zero = FALSE,
                                trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16)
  p.format <- stringr::str_to_sentence(p.format)
  
  # associated model evaluation metrics 
  lab <- paste(paste0("MAE=", round(evalutaion$MAE,3)), 
               paste0("R=", round(evalutaion$Cor,3)),
               p.format, sep="; ")
  
  # data for slope line
  dat.mm <- data.frame(min=min(prediction$obs), max=max(prediction$obs))
  
  p1 <- ggplot(prediction, aes(obs, pred)) +
    geom_segment(aes(x=min, y=min, xend=max, yend=max),
                 lty="dashed", data=dat.mm) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic() + 
    theme(plot.title = element_text(hjust=0.5, face = "bold"),
          plot.subtitle = element_text(hjust=0.5),
          axis.text = element_text(color="black")) +
    labs(x="Chronological Age", y="Predicted Age", 
         subtitle=lab, title=title)
  
  return(p1)
}

# predict with validation set
validPred <- function(model, obs,
                      nad.mat=NULL, expr.mat=NULL, 
                      model.coef=NULL
) {
  
  # Model coefficients  
  if (is.null(model.coef)) {
    coef_model <- rownames(coef(model))[-1]  
  } else {
    coef_model <- model.coef
  }
  
  
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

# prediction with leave-one-out (LOO) 
LOOP <- function(data, y, model=c("glmnet","lasso","ridge","randomforest"),
                 alpha=NULL, lambda=NULL, mtry=NULL) {
  # use all data for training
  pred.var <- model.matrix(formula(paste(y, "~ .")), data=data)[,-1]
  resp.var <- data[, y]
  
  # Create an empty vector to store the LOOCV predictions
  cv.pred <- rep(0,length(resp.var))
  idx <- 1:length(resp.var)
  # Perform LOO analysis
  for (j in idx) {
    # Split the data into training and test sets
    train.index <- idx[-j]
    test.index <- idx[j]
    train.data <- pred.var[train.index,]
    train.labels <- resp.var[train.index]
    test.data <- pred.var[test.index,]
    
    if (model %in% c("glmnet", "lasso", "ridge")) {
      fit.model <- glmnet(train.data, train.labels, 
                          family="gaussian",
                          alpha = alpha, 
                          lambda = lambda)
    } else if (model == "randomforest") {
      # Fit the random forest model using the training set
      fit.model <- randomForest(train.data, train.labels, ntree=1000, mtry=mtry)
    }
    
    # Make predictions using the test set
    cv.pred[j] <- predict(fit.model, test.data)
  }
  
  pred.df <- data.frame(obs=resp.var, pred=cv.pred)
  eval.mat <- evalMetrics(pred.df$obs, pred.df$pred)
  
  return(list(prediction=pred.df, evaluation=eval.mat))
}

# use when: “Error in summary.connection(connection) : invalid connection”
# https://stackoverflow.com/questions/64519640/error-in-summary-connectionconnection-invalid-connection
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# dual barplot
barplot_pubr <- function(enrich.obj, 
                         term_col="Description", 
                         term_size_col="Count", 
                         p_col="p.adjust",
                         showCategory=10, title=NULL, str_width=35, 
                         log=TRUE){
  eobj <- enrich.obj[, c(term_col, term_size_col, p_col)]
  
  eggtmp <- head(eobj, n=showCategory)
  
  if (log) {
    eggtmp[, p_col] <- -log10(eggtmp[, p_col])
  }
  
  eggtmp <- reshape2::melt(eggtmp)
  eggtmp[, term_col] <- factor(eggtmp[, term_col], levels=rev(unique(eggtmp[, term_col])))
  
  bp <- ggplot(eggtmp, aes_string(term_col)) +
    geom_bar(aes(y=value, fill=variable), stat='identity', 
             position=position_dodge(width=0.9), width=0.8) +
    coord_flip() +
    theme_minimal() +
    geom_hline(yintercept = (-log10(0.05)),lty=4,col="grey",lwd=0.6) +
    theme(axis.text.y = element_text(color='black'),
          axis.text.x = element_text(color='black'),
          axis.line = element_line(),
          axis.ticks = element_line(color='black'),
          legend.position = 'right',
          panel.grid = element_blank()) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = str_width)) +
    scale_fill_manual(values = c("#2d74b9","#d28b47"), labels = c('Gene Counts', '-Log10(p.adjust)')) +
    labs(x='',y='',fill='', title=title)
  return(bp)
}

# convert p value to significance mark
# https://stackoverflow.com/questions/41262992/is-there-a-r-function-that-convert-p-value-to-significance-code
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
         symbols = c("***", "**", "*", ""))
}

# rank metrics
rankMetrics <- function(metrics) {
  # multiplying by +/- 1 so that large values correspond to good performance
  score.mat <- t(t(metrics) * c(1,1,-1,1,-1,-1,1,-1)) # BIO_SIL,ASSAY_SIL,BATCH_SIL,PAM_SIL,RLE_MED,RLE_IQR,EXP_WV_COR,EXP_UV_COR
  rank.mat <- apply(na.omit(score.mat), 2, rank, ties.method='min')
  # if score all 1, remove it before scoring
  metrics.keep <- colSums(rank.mat==1) != nrow(rank.mat)
  mean.score <- rowMeans(rank.mat[,metrics.keep])
  rank.mat <- as.data.frame(rank.mat)
  rank.mat$SCORE <- mean.score
  return(rank.mat)
}

# synthetic RNA helper
synHelp <- function(object, method="TC", log=TRUE) {
  
  if (!any(SummarizedExperiment::rowData(object)$Synthetic)) {
    stop("Synthetic RNA id not provided in the object")
  }
  
  # get raw counts
  counts_df <- SummarizedExperiment::assay(object)
  
  # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
  method.curr <- unlist(strsplit(method, split = "_"))
  
  # check whether selected method can be provided
  if (!method.curr[1] %in% c("Raw","TC","UQ","DESeq","TMM","PossionSeq")) {
    stop("Scaling method: ", method.curr[1], " is not provided. It should be one of ", c("Raw","TC","UQ","DESeq","TMM","PossionSeq"))
  }
  if (length(method.curr) > 1 & !method.curr[2] %in% c("RUVg","RUVs","RUVse")) {
    stop("RUV method: ", method.curr[2], " is not provided. It should be one of ", c("RUVg","RUVs","RUVse"))
  }
  if (length(method.curr) > 2 & is.na(stringr::str_extract(method.curr[3],"k"))) {
    stop("Use x number of factors to estimate unwanted variation as \"kx\", but not ", method.curr[3])
  }  
  if (length(method.curr) > 2 & as.numeric(stringr::str_extract(method.curr[3],"\\d")) > ncol(object)) {
    stop("Number of required factors exceed, try least.")
  } 
  
  # scaling
  if (method.curr[1] == "Raw") {
    counts_scale <- list(dataNorm=counts_df, normFactor=rep(1,ncol(counts_df)))
  } else {
    normScaling <- get(paste0("norm", method.curr[1]))
    counts_scale <- normScaling(counts_df)
  }
  
  # RUV
  sc_idx <- CreateGroupMatrix(object$condition)
  enrich_idx <- CreateGroupMatrix(object$enrich)
  if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
    
    sc.idx <- switch(method.curr[2],
                     "RUVg" = NULL,
                     "RUVs" = sc_idx,
                     "RUVse" = enrich_idx)
    counts_norm <- normRUV(counts_scale$dataNorm,
                           control.idx = getGeneSet(object, "NegControl"),
                           sc.idx = sc.idx,
                           method = method.curr[2],
                           k = as.numeric(gsub("k", "", method.curr[3])))
    
  } else {
    counts_norm <- counts_scale
  }
  
  # calculate enrichment of synthetic RNA
  counts_norm <- counts_norm$dataNorm
  syn_id <- SummarizedExperiment::rowData(object)$Synthetic
  # add 1 offset to avoid zero division (in log-scale)
  syn_en <- log2(counts_norm[syn_id, enrich_idx[1,]]+1) - log2(counts_norm[syn_id, enrich_idx[2,]]+1)
  
  if (log) {
    syn_en <- syn_en
  } else {
    syn_en <- 2^syn_en
  }
  
  return(list(enrichment=syn_en, dataNorm=counts_norm[syn_id,]))
}

# simple RLE (relative log-expression) transformer
expr2rle <- function(x) {
  x.log <- log2(x+1)
  x.log.rle <- x.log - median(x.log)
}

# normalization from RADAR
normRADAR <- function(data, en.idx, top.enrich = 0.01) {
  # get input data
  dataInput <- data[, en.idx[2,]]
  # estimate size factors for input
  sizeFactorInput <- DESeq2::estimateSizeFactorsForMatrix(dataInput)
  # normalized input
  InputNorm <- t(t(dataInput)/sizeFactorInput)
  
  # get top k% enrich genes
  dataEnrich <- data[, en.idx[1,]]
  avgEnrich <- sort(rowMeans(dataEnrich), decreasing = TRUE)
  topEnrich.id <- names(head(avgEnrich, n=length(avgEnrich)*top.enrich))
  
  # calculate enrichment based on top enrichment
  topEnrich <- dataEnrich[topEnrich.id,]/InputNorm[topEnrich.id,]
  
  # estimate size factors for enrich
  sizeFactorEnrich <- DESeq2::estimateSizeFactorsForMatrix(topEnrich)
  
  # normalized Enrich
  EnrichNorm <- t(t(dataEnrich)/sizeFactorEnrich)
  
  # compute gene-wise size factor to adjust pre-IP gene expression
  geneSF <- InputNorm/rowMeans(InputNorm)
  geneSF[geneSF==0] <- 1
  # adjust Enrich
  EnrichAdjust <- EnrichNorm/geneSF
  
  dataNorm <- cbind(InputNorm, EnrichAdjust)[,colnames(data)]
  return(dataNorm)
}

# ARI
calcARI <- function(label, cluster) {
  pair_mat <- table(label, cluster)
  
  # RI: number of success
  ri <- sum(apply(pair_mat, 2, function(x) {
    # for each succeed pairs
    col.k <- sum(sapply(x[x>1], choose, k=2))
  }))
  # for each row
  a <- sum(sapply(rowSums(pair_mat), choose, k=2))
  # for each column
  b <- sum(sapply(colSums(pair_mat), choose, k=2))
  # number of pairs
  n_pair <- choose(length(label), 2)
  
  # eRI
  eri <- a*b/n_pair
  # maxRI
  maxri <- (a+b)/2
  
  ARI <- (ri - eri)/(maxri - eri)
  return(ARI)
}

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
    # annotate(geom="text", x = 25, y = -1.4, 
    #          label=paste0("n = ", length(unique(unlist(data.plot[, group]))))) +
    theme_classic() +
    theme(axis.text = element_text(color="black")) +
    labs(x="Age", y="Z-score", title=title)
  
  gg1 
}

##----##