## Function to retrieve Mscore for one patient and one path against healthy control distribution
.getMscorePath <- function(path,
                         Patient,
                         Healthy){
  genes.path <- Reduce(intersect, list(rownames(Healthy),
                                    names(Patient),
                                    as.character(unlist(path))))
  path.mscore <- 0 # default

  if(length(genes.path) > 2){ # at least 3 genes in each path
    tmpRef <- Healthy[genes.path,]
    tmpPat <- Patient[genes.path]

    # Zscore by path
    sel <- ifelse(apply(tmpRef, 1, sd) == 0, FALSE, TRUE)

    if(sum(sel) > 2){ # only in modules with more than two genes
        tmpRef <- tmpRef[sel,]
        tmpPat <- tmpPat[sel]

        Zscore.genes <- (tmpPat-apply(tmpRef, 1, mean))/apply(tmpRef, 1, sd)
        path.mscore <- mean(Zscore.genes, na.rm=TRUE)
    }
  }
  names(path.mscore) <- names(path)
  return(path.mscore)
}







#-------------------------------------------------------------- .getNearSample()
## Impute m-score for a patient from gene expression based on a gene expression reference of patients
#@ patient: vector of gene expression of a patient (names of vector must be gene symbols)
#@ Reference.normalized: gene expression matrix of patients normalized by z-score (by patient)
#@ Reference.mscore: matrix with mscores of patients contained in @Reference.normalized
#@ k: number of neighbours to impute mscore for @patient
.getNearSample <- function(patient,
                         Ref.norm,
                         Ref.mscore,
                         k=5){
  patient <- .normSamples(patient[intersect(names(patient), rownames(Ref.norm))])

  distances <- apply(Ref.norm, 2, function(x) {
    dist(rbind(patient, x))
  })
  names(distances) <- colnames(Ref.norm)

  distances <- distances[order(distances,decreasing=FALSE)]
  tmp.mscore <- Ref.mscore[,names(distances)[1:k]]
  tmp.mscore <- apply(tmp.mscore,1,mean)

  return(tmp.mscore)
}


#----------------------------------------------------------------- .normSamples()
## Normalize quantitative values (expression, mscores) of a patient by z-score
#@ x: Numeric vector of expression from a sample
.normSamples <- function(x){
  x <- (x-(mean(x,na.rm=TRUE)))/sd(x,na.rm=TRUE)
  return(x)
}


#------------------------------------------------------------------- .fixSymbol()
## Modify character patterns within a character vector
#@ vect: character vector
#@ symbol: patterns to replace
#@ toSymbol: patterns with which the patterns of @symbol are replaced (position by position)
.fixSymbol <- function(vect,
                    symbol=c("/"," "),
                    toSymbol=c(".","")){
  for(i in 1:length(symbol)){
    vect <- gsub(pattern=symbol[i],replacement=toSymbol[i],vect)
  }
  return(vect)
}


#------------------------------------------------------------------ .MReference()
## Function to build a normalized reference of patients with which to impute the m-scores for new samples
#@ expr.list: list with expression data from one or multiple studies (used as reference)
#@ geneset.list: output of DiseasePaths function (Pathways relevant for the disease)
#@ mscore.list: list with m-score matrices from the same studies that expr.list
.MReference <- function(expr.list,
                     geneset.list,
                     mscore.list){

  ## Get common genes for all samples of the reference
  all.genes <- lapply(expr.list, function(i){x <- rownames(i[[1]])})
  common.genes <- Reduce(intersect, all.genes)

  # COMPROBAR MÉTODO DANI
  # all.genes[["genesets"]] <- unique(as.character(unlist(geneset.list)))
  # common.genes <- all.genes[[1]]
  # for (i in 2:length(all.genes)) {
  #   common.genes <- intersect(as.character(all.genes[[i]]), common.genes)
  # }

  ## Reference gene-expression
  Reference <- lapply(expr.list, function(x) {return(x[[1]][common.genes,])})
  Reference <- do.call(cbind, Reference)
  Reference.normalized<-apply(Reference, 2, .normSamples)

  ## Reference M-score
  Reference.mscore <- lapply(mscore.list, function(x) {
      return(x[names(geneset.list),])})
  Reference.mscore <- do.call(cbind, Reference.mscore)

  return(list(Reference.mscore=Reference.mscore,
              Reference.normalized=Reference.normalized))
}


#----------------------------------------------------------------------- GetML()
## Build ML-based model to predict clinical outcomes from expression/mscores
#@ exp.data: expression/mscore matrix. Samples in columns and features in rows
#@ metadata: dataframe with information for each sample. Samples in rows and variables in columns
#@ var2predict: character with the column name of the @metadata to predict
#- Values of metadata$var2predict must be "YES"/"NO" (for categoral outcomes) or quantitative (for numerical outcomes)
#@ var.type: type of outcome to predict (categorical (cat) or numerical (num))
#@ algorithms: algoritms for ML (glm,lm,lda,xgbTree,rf,knn,svmLinear,svmRadial,nnet,nb,lars,rpart). "all" = all algorithms are used
#@ outerfolds: number of external folds in which exp.data is divided
#@ prop: proportion of samples used as train for the outerfolds (0-1)
#@ innerfold: number of internal folds (for parameter tuning)
#@ repeats: number of repetitions of the parameter tuning process
#@ feature.filter: method to reduce number of features (none,fcbf)
#@ prior: rank best model based on AUC (AUC), the mean of AUC and the balanced accuracy (BAUC), the balanced accuracy (BA) or manual (Manual)
## FALTA PONER PARÁMETRO ADD PARA AÑADIR MODELOS EXTRA
GetML <- function(exp.data,
                metadata,
                var2predict,
                var.type="cat",
                algorithms="all",
                outerfolds=10,
                prop=0.8,
                innerfold=10,
                repeats=10,
                feature.filter="none",
                prior="AUC"){
  check.packages(c("caret","mlbench","pROC","rpart","randomForest","nnet",
                   "caretEnsemble","MLeval","pROC","ROCR","ada", "plyr",
                   "mboost","import","brnn","elasticnet"))

  if(var2predict %in% colnames(metadata)){

    samples <- intersect(rownames(metadata),colnames(exp.data))
    exp.data <- exp.data[,samples]
    metadata <- metadata[samples,]

    train.data <- data.frame("group"=metadata[,var2predict],
                           as.data.frame(t(exp.data)))
    colnames(train.data) <- .fixSymbol(vect=colnames(train.data),
                                    symbol=c("/"," ","-"),
                                    toSymbol=c(".",".","."))
    train.data <- train.data[ifelse(is.na(train.data$group),F,T),]
    #var.type=ifelse(is.character(train.data$group[1]),"cat","num")

    result.ML <- .neastkfoldML(data=train.data,
                            var.type=var.type,
                            algorithms=algorithms,
                            outerfolds=outerfolds,
                            kfold=innerfold,
                            prop=prop,
                            repeats=repeats,
                            feature.filter=feature.filter,
                            prior=prior)

    return(result.ML)
  }else{
    error <- "Variable to predict not found in metadata"
    return(error)
  }
}


#-------------------------------------------------------------------- .MLmethod()
## Prepare parameters for tuning the different methods
#@ algorithms: algoritms for ML (glm,lm,lda,xgbTree,rf,knn,svmLinear,svmRadial,nnet,nb,lars,rpart). "all" = all algorithms are used
#@ var.type: type of outcome to predict (categorical (cat) or numerical (num))
#@ training: train data
#@ add: additional model list (supported by caret) in format: add=list(modelName = caretModelSpec(method="modelName", tuneGrid=(.parameters="values")))
.MLmethod <- function(algorithms,
                   var.type,
                   training,
                   add=NULL){
  if("all" %in% algorithms){
    algorithms <- c("glm", "lm", "lda", "xgbTree", "rf", "knn", "svmLinear", "svmRadial",
                  "nnet", "nb", "lars", "rpart", "ada", "gamboost", "brnn", "enet")}

  switch(var.type,
         cat={
           algorithms <- algorithms[algorithms %in% c("glm", "lda", "xgbTree", "rf", "knn", "svmLinear", "svmRadial", "nnet", "nb", "ada", "gamboost")]
         },
         num={
           algorithms <- algorithms[algorithms %in% c("glm", "lm", "xgbTree", "rf", "knn", "svmLinear", "svmRadial", "nnet", "lars", "rpart", "gamboost", "brnn", "enet")]
         })

  method.list <- list(
    lm=caretModelSpec(method="lm"),
    glm=caretModelSpec(method="glm"),
    lda=caretModelSpec(method="lda"),
    xgbTree=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(.max_depth=c(2, 3, 4, 5, 6, 8, 10), .nrounds=50, .eta=c(0.01, 0.05, 0.1), .gamma=c(0, 1), .colsample_bytree=c(0.1,  0.4),  .min_child_weight=c(1,  10),  .subsample=c(0.5,  1))),
    rf=caretModelSpec(method="rf", tuneGrid=data.frame(.mtry=seq(2, round(sqrt(ncol(training)-1)), 1))),
    knn=caretModelSpec(method="knn", tuneGrid=data.frame(.k=seq(2, 50, 2))),
    svmLinear=caretModelSpec(method="svmLinear", tuneLength=15),
    svmRadial=caretModelSpec(method="svmRadial", tuneLength=15),
    nnet=caretModelSpec(method="nnet", tuneGrid=expand.grid(.size=seq(from=1, to=4, by=1),.decay=seq(from=0.1, to=0.5, by=0.1))),
    nb=caretModelSpec(method="nb", tuneGrid=expand.grid(.fL=c(0, 0.5, 1.0), .usekernel=TRUE, .adjust=c(0.5, 1.0))),
    lars=caretModelSpec(method="lars", tuneGrid=expand.grid(.fraction=seq(.01, .99, length=40))),
    rpart=caretModelSpec(method="rpart", tuneGrid=expand.grid(.cp=seq(0, .02, .0001))),
    ada=caretModelSpec(method="ada", tuneGrid=expand.grid(.maxdepth=25, .nu=2, .iter=100)),
    gamboost=caretModelSpec(method="gamboost", tuneGrid=expand.grid(.mstop=c(50, 100, 150, 200, 250, 300), .prune=c('yes', 'no'))),
    brnn=caretModelSpec(method="brnn", tuneGrid=expand.grid(.neurons=c(2:16))),
    enet=caretModelSpec(method="enet", tuneGrid=expand.grid(.lambda=c(0, 0.005, 0.009, 0.05, 0.1), .fractions=c(0.05, 0.09, 0.5)))
  )

  method.list <- method.list[algorithms]
  method.list <- c(method.list,add)

  return(method.list) # methodList and tuneList
}


#---------------------------------------------------------------- .neastkfoldML()
## ## Get performance results for different algorithms
#@ data: matrix with features in columns and samples in rows. Column of group is required (outcome to predict)
#@ algorithms: algoritms for ML (glm,lm,lda,xgbTree,rf,knn,svmLinear,svmRadial,nnet,nb,lars,rpart). "all" = all algorithms are used
#@ outerfolds: number of external folds in which exp.data is divided
#@ kfold: innerfold; number of internal folds (for parameter tuning)
#@ repeats: number of repetitions of the parameter tuning process
#@ prop: proportion of samples used as train for the outerfolds (0-1)
#@ var.type: type of outcome to predict (categorical (cat) or numerical (num))
#@ feature.filter: method to reduce number of features (none,fcbf)
#@ prior: rank best model based on AUC (AUC), the mean of AUC and the balanced accuracy (BAUC), the balanced accuracy (BA) or manual (Manual)
.neastkfoldML <- function(data,
                       algorithms,
                       outerfolds,
                       kfold,
                       repeats,
                       prop,
                       var.type,
                       feature.filter,
                       prior="BAUC"){
  check.packages(c("caret", "mlbench", "pROC", "rpart", "randomForest", "nnet", "caretEnsemble", "MLeval", "pROC", "ROCR"))

  ##-------------- (STEP 1)
  ## Outer folds

  list.perm.bias <- list()
  for(pb in 1:outerfolds){ ## Create outer folds
    if(outerfolds<2){outerfolds=2}
    list.perm.bias[[pb]] <- createDataPartition(y=data$group, p=prop, list=FALSE)
  }

  ## 1.1 Run outer fold
  outerML <- list()
  for(i in 1:length(list.perm.bias)){

    ## Train/ test
    training <- data[list.perm.bias[[i]],]
    testing <- data[-list.perm.bias[[i]],]

    my_control <- trainControl(
      method="repeatedcv",
      number=kfold, ## inner folds
      savePredictions="final",
      classProbs=TRUE,
      index=createResample(training$group, kfold),
      repeats=repeats)

    ## 1.2 Feature selection
    switch(feature.filter,
           "fcbf"={
             cat("\nSelecting features by fcbf")
             training <- .fastCorFS(training, thresh=0.0025)
             testing <- testing[, colnames(training)]
           },
           "none"={
             cat("\nNo feature selection selected")
           })

    ## 1.3 Parameter tuning
    if(var.type=="cat"){ ## Categorical outcome

      ## Get ML-methods
      method.list <- .MLmethod(algorithms=algorithms, var.type=var.type, training=training)

      model_list <- caretList(
        group~., data=training,
        trControl=my_control, ## inner folds
        #methodList=c("glm","lda"),
        tuneList=method.list)

      ## Get prediction and performance in outer test (for algorithm selection)
      model_pred <- list()
      cm <- list()
      for(mp in 1:length(model_list)){

        p <- predict(model_list[[mp]], newdata=testing, type="prob")
        if(class(p)=="numeric"){
          cm.model <- confusionMatrix(as.factor(ifelse(p<=0.5, "YES", "NO")), as.factor(testing$group))
          pred <- data.frame("NO"=p, "YES"=1-p, "obs"=testing$group)
        }else{ ## data.frame
          cm.model <- confusionMatrix(as.factor(ifelse(p$NO<=0.5, "YES", "NO")), as.factor(testing$group))
          pred <- data.frame("NO"=p$NO, "YES"=1-p$NO, "obs"=testing$group)
        }
        pred <- na.omit(pred)
        model_pred[[mp]] <- pred
        cm[[mp]] <- cm.model
      }
      names(model_pred) <- names(model_list)
      res <- list(model_list, model_pred, cm)
      names(res) <- c("models", "preds", "cm")

      outerML[[i]] <- res

    }else { ## numerical outcome

      ## Get ML-methods
      method.list <- .MLmethod(algorithms=algorithms, var.type=var.type, training=training)

      model_list <- caretList(
        group~., data=training,
        trControl=my_control, ## inner folds
        #methodList=c("glm","lm"),
        tuneList=method.list)

      ## Get prediction and performance in outer test (for algorithm selection)
      model_pred <- list()
      for(mp in 1:length(model_list)){
        p <- predict(model_list[[mp]], newdata=testing)
        model_pred[[mp]] <- data.frame("pred"=as.numeric(p), "obs"=as.numeric(testing$group))
      }
      names(model_pred) <- names(model_list)

      res <- list(model_list, model_pred)
      names(res) <- c("models", "preds")
      outerML[[i]] <- res
    }
  } # outer folds

  ##-------------- (STEP 2)
  ## Best algorithm selection

  ## 2.1 Get performance metrics for all algorithms
  if(var.type=="cat"){ ## Categorical outcome
    STATS <- matrix(data=0, ncol=length(outerML[[1]]$models), nrow=13)
    colnames(STATS) <- names(outerML[[1]]$models)

    for(mod in 1:ncol(STATS)){

      tmp <- outerML[[1]]$preds[[mod]]
      cm.model <- outerML[[1]]$cm[[mod]]
      stats <- c(cm.model$overall["Accuracy"],
               cm.model$byClass[c("Balanced Accuracy", "Precision", "Recall", "F1",
                                  "Prevalence", "Detection Rate", "Detection Prevalence",
                                  "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")])

      count <- ifelse(is.na(stats), 0, 1)
      for(i in 2:length(outerML)){
        tmp <- rbind(tmp, outerML[[i]]$preds[[mod]])

        cm.model <- outerML[[i]]$cm[[mod]]
        stats2 <- c(cm.model$overall["Accuracy"],
                  cm.model$byClass[c("Balanced Accuracy", "Precision", "Recall", "F1",
                                     "Prevalence", "Detection Rate", "Detection Prevalence",
                                     "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")])

        count <- count+ifelse(is.na(stats2), 0, 1)
        stats2[is.na(stats2)] <- 0
        stats <- stats+stats2
      }
      stats <- stats/count

      x <- NA
      tryCatch({
        x <- evalm(tmp)
        x <- x[[7]][[1]]["AUC-ROC","Score"]
      }, error=function(e){})
      stats <- c(x, stats)
      names(stats)[1] <- "AUC"
      STATS[, mod] <- stats
    }
    rownames(STATS) <- names(stats)

    ## 2.2 Select best algorithm
    if (ncol(STATS) > 1) {
      switch(prior,
             "AUC"={
               STATS <- STATS[, order(STATS["AUC",], decreasing=TRUE)]
             },
             "BAUC"={
               ord <- apply(STATS[c(1, 3), ], 2, sum)
               STATS <- STATS[, order(ord, decreasing=TRUE)]
             },
             "Manual"={
               print(STATS)
               sel <- readline("Select algorithm > ")
               sel <- as.character(sel)
               STATS <- STATS[, c(sel, setdiff(colnames(STATS), sel))]
             },
             "BA"={
               STATS <- STATS[, order(STATS["Balanced Accuracy", ], decreasing=TRUE)]
             })
    }


    ## 2.3 Build model with all data, best algorithm and best parameters
    bestTune <- outerML[[1]]$models[[colnames(STATS)[1]]]$bestTune
    for(i in 2:length(outerML)){
      bestTune <- rbind(bestTune, outerML[[i]]$models[[colnames(STATS)[1]]]$bestTune)
    }
    bestTune <- apply(bestTune, 2, mean)
    bestTune <- as.data.frame(t(bestTune), nrow=1)
    colnames(bestTune) <- paste0(".", colnames(bestTune))

    ## Feature selection
    switch(feature.filter,
           "fcbf"={
             data.f <- .fastCorFS(data, thresh=0.0025)
           },
           "none"={
             data.f <- data
           })

    fit.model <- train(group ~ ., data=data.f, method=colnames(STATS)[1],
                     tuneGrid=bestTune)
    stats <- STATS
    preds <- NULL
    for(pr in 1:length(outerML)){
      preds <- rbind(preds, outerML[[pr]]$preds[[colnames(STATS)[1]]])
    }
    auc <- evalm(preds)

    model <- list(fit.model, stats, auc, bestTune)
    names(model) <- c("model", "stats", "auc", "bestTune")

    return(model)

  }else{ ## Numerical outcome

    ## 2.1 Get performance metrics for all algorithms
    STATS <- rep(0, length(outerML[[1]]$models))
    names(STATS) <- names(outerML[[1]]$models)

    for(mod in 1:length(outerML[[1]]$models)){
      tmp <- outerML[[1]]$preds[[mod]]
      for(pr in 2:length(outerML)){
        tmp <- rbind(tmp, outerML[[pr]]$preds[[mod]])
      }
      STATS[mod] <- cor(tmp[, 1], tmp[, 2])
    }
    ## 2.2 Select best model by correlation
    STATS <- STATS[order(STATS, decreasing=TRUE)]

    ## 2.3 Build model with all data, best algorithm and best parameters
    bestTune <- outerML[[1]]$models[[names(STATS)[1]]]$bestTune
    for(i in 2:length(outerML)){
      bestTune <- rbind(bestTune, outerML[[i]]$models[[names(STATS)[1]]]$bestTune)
    }

    bestTune <- apply(bestTune, 2, mean)
    bestTune <- as.data.frame(t(bestTune), nrow=1)
    colnames(bestTune) <- paste0(".", colnames(bestTune))

    ## Feature selection
    switch(feature.filter,
           "fcbf"={
             data.f <- .fastCorFS(data, thresh=0.0025)
           },
           "none"={
             data.f <- data
           })

    fit.model <- train(group ~ ., data=data.f, method=names(STATS)[1],
                     tuneGrid=bestTune)
    stats <- STATS

    model <- list(fit.model, stats, bestTune)
    names(model) <- c("model", "stats", "bestTune")

    return(model)
  }
}


#----------------------------------------------------------------- .fastCorFS()
## Function to remove co-linear features
#@ data: matrix with features in columns and samples in rows. Column of group is required (outcome to predict)
#@ thresh: SU threshold
.fastCorFS <- function(data,
                      thresh=0.0025){
  stopifnot('group' %in% colnames(data))
  check.packages("FCBF")

  y <- as.factor(data$group)
  x <- subset(data, select=-c(group))
  dis <- discretize_exprs(t(x))

  fcbf.res=fcbf(dis, y, verbose=TRUE, thresh)
  xx <- x[, fcbf.res$index]
  xx=as.data.frame(cbind(group=y, xx))

  return(xx)
}


#---------------------------------------------------------------------- cuteML()
## Minimun-Feature selection and normalize ML models (only for models that support varImp function)
#@ listfit: model, output object of GetML function
#@ norm: normalize values using z-score by patients (recommended for gene-expression based models, not if m-scores are used)
#@ sel: Character indicates how to select best model (max: number of features that maximize the performance, manual: manual)
#@ nmax: maximun of pathways to test
#@ pathPlot: path to save output plots
#@ seed: set.seed, number
#@ ...: additional parameters for the plot (res=300, width=3, height=1.5, units="in")
## QUITAR NORM
cuteML <- function(listfit,
                 norm=FALSE,
                 sel="max",
                 nmax=NULL,
                 pathPlot=NULL,
                 seed=12345678,
                 ...){
  check.packages(c("scales", "caret", "ggplot2"))
  set.seed(seed)

  data <- listfit$model$trainingData ## extract data from the model object
  colnames(data)[1] <- "group"

  ## Feature Importance FALTA (buscar alternativa a varImp para el resto de modelos y numericals)
  imp <- cbind(varImp(listfit$model)$importance,
             varImp(listfit$model)$importance)
  imp <- imp[order(imp[, 1], decreasing=TRUE),]

  nmax <- ifelse(is.null(nmax), nrow(imp), nmax)
  nmin <- ifelse(norm, 2, 1)
  acc <- NULL
  models <- list()
  count <- 1
  for(i in nmin:nmax){
    genes <- c("group", rownames(imp)[1:i])
    tmp <- data[, genes]

    if(norm){
      tmp.z <- tmp
      for(z in 1:nrow(tmp.z)){
        tmp.z[z, 2:ncol(tmp.z)] <- as.numeric(scale(as.numeric(tmp[z, 2:ncol(tmp)]), center=TRUE, scale=TRUE))
      }
      tmp <- tmp.z
    }

    mod.i <- train(group ~ ., data=tmp, method=listfit$model$method,
                 #trControl=my_control,
                 tuneGrid=listfit$bestTune)
    models[[count]] <- mod.i
    names(models)[length(models)] <- i
    count <- count+1
    acc <- c(acc, mod.i$results$Accuracy)
  }

  ## Plot
  m <- data.frame("nfeatures"=nfeatures <- nmin:nmax,
                "accuracy"=acc)
  p1 <- ggplot(m, aes(x=nfeatures, y=accuracy))+theme_light()+
    geom_point(size=0.5)+geom_line()+
    theme(axis.text.x=element_text(size =5, angle=90, vjust=0.5, hjust=1),
          axis.text.y=element_text(size=5),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    scale_x_continuous(breaks=round(seq(min(m$nfeatures), max(m$nfeatures), by=2), 1))
  if(!is.null(pathPlot)){
    tiff(filename=pathPlot, ...)
    plot(p1)
    invisible(dev.off())
  }

  ## Minimun-Feature selection
  switch(sel,
         "max"={
           nfeatures <- m[which.max(m$accuracy), "nfeatures"]
         },
         "manual"={
           plot(p1)
           nfeatures <- as.numeric(readline("Number of features > "))
         })

  model.final <- models[as.character(nfeatures)][[1]]

  return(model.final)
}
