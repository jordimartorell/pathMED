#' Build machine-learning based on gene-expression/mscores data
#'
#' @param expData Feature matrix. Samples in columns and features in
#' rows. expData must be a numerical matrix
#' @param metadata dataframe with information for each sample. Samples in rows
#' and variables in columns
#' @param models named list with the ML models generated with
#' caret::caretModelSpec function. methodsML function may be used to prepare
#' this list
#' @param var2predict character with the column name of the @metadata to predict
#' @param paired character with the column name of the @metadata with paired information.
#' Paired information is used to always asign samples from the same source (i.e. patient)
#' to the same group (train or test). paired = NULL by default
#' @param Koutter number of Outter folds in which exp.data is split for neasted
#' cold validation.
#' A list of integer with elements for each resampling iteration is admitted.
#' Each list element is a vector of integers corresponding to the rows used for
#' training at that iteration.
#' @param Kinner number of cross-validation folds (for parameter tuning)
#' @param repeatsCV number of repetitions of the parameter tuning process
#' @param filterFeatures "rfe" (Recursive Feature Elimination), "sbf" (Selection
#' By Filtering) or NULL (no feature selection)
#' @param filterSizes Only for filterFeatures = "rfe". A numeric vector of integers
#' corresponding to the number of features that should be retained
#' @param rerank Only for filterFeatures = "rfe". A logical: should variable importance be
#' re-calculated each time features are removed?
#' @param positiveClass outcome value that must be considered as positive class
#' (for categoric outcomes)
#' @param continue_on_fail whether or not to continue training the models if any of them fail
#' @param saveLogFile path to a .txt file in which to save error and warning messages
#'
#' @return A list with four elements. The first one is the model. The second one
#' is a table with different metrics obtained. The third one is a list with the
#' best parameters selected in tuning process. The last element contains data
#' for AUC plots
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @import caret
#' @import caretEnsemble
#' @import xgboost
#' @import randomForest
#' @import klaR
#' @import ada
#' @import mboost
#' @import import
#' @import kernlab
#' @import metrica
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(refData, exampleData, exampleMetadata)
#' \donttest{
#' exampleRefMScore <- createReference(data=refData, genesets="tmod")
#' relevantPaths <- diseasePaths(MRef=exampleRefMScore,
#' min_datasets=3,
#' perc_samples=10)
#'
#' MScoresExample <- getMscores(genesets = relevantPaths,
#' Patient = exampleData,
#' Healthy = NULL,
#' nk = 5)
#'
#' fit.model <- getML(expData=MScoresExample,
#' metadata=exampleMetadata,
#' var2predict="Response",
#' models=methodsML("svmLinear"),
#' Koutter=5,
#' Kinner=10,
#' repeatsCV=10)
#' }
#'
#' @export
getML <- function(expData,
                  metadata,
                  models=methodsML(outcomeClass = "character"),
                  var2predict,
                  paired=NULL,
                  Koutter=5,
                  Kinner=4,
                  repeatsCV=5,
                  filterFeatures=NULL,
                  filterSizes = seq(2,100, by = 2),
                  rerank=FALSE,
                  continue_on_fail = TRUE,
                  positiveClass=NULL,
                  saveLogFile = NULL){
  if (!is.null(saveLogFile)) {
    if (file.exists(saveLogFile)) {
      file.remove(saveLogFile) # reset log file
    }
  }

  if (!var2predict %in% colnames(metadata)) {
    stop("var2predict must be a column of metadata")
  }

  if(!is.null(filterFeatures)){
    if(!filterFeatures %in% c("sbf","rfe")){
      stop("filterFeatures must be 'sbf', 'rfe', or NULL")
    }else{
      if(filterFeatures=="rfe"){
        filterSizes <- filterSizes[filterSizes <= nrow(expData)]
      }
    }
  }
  ## 1. Formatting data input
  samples <- intersect(rownames(metadata), colnames(expData))

  if (length(samples) < 1) {
    stop("Row names of metadata and column names of expData do not match")
  }

  if (is.character(metadata[, var2predict])|is.factor(metadata[, var2predict])) {
    metadata[, var2predict] <- stringi::stri_replace_all_regex(metadata[, var2predict],
                                                               pattern = c("/", " ", "-"), replacement = c(".", ".","."), vectorize = FALSE)
  }
  if (!is.null(positiveClass)) {
    positiveClass <- stringi::stri_replace_all_regex(positiveClass,
                                                     pattern = c("/", " ", "-"), replacement = c(".", ".","."), vectorize = FALSE)
  }

  expData <- expData[,samples]
  metadata <- metadata[samples,,drop=FALSE]

  expData <- expData[rowSums(expData) != 0, ] # Remove features with all 0

  expData <- data.frame('group'=metadata[,var2predict],
                        as.data.frame(t(expData)))
  colnames(expData) <- stringi::stri_replace_all_regex(
    colnames(expData), pattern = c('/',' ','-',':'),
    replacement = c('.','.','.','.'), vectorize=FALSE)
  expData <- expData[!is.na(expData$group),]
  if(is.factor(expData$group)){
    expData$group<-as.character(expData$group)
  }
  if(expData %>% summarise(across(everything(), ~ any(is.na(.) | is.infinite(.)))) %>% any()){
    stop("There are NA and/or Infinite values in your data, please replace or remove them before running getML")
  }

  outcomeClass <- class(expData$group)
  if (methods::is(expData$group, "character")){
    prior <- "MCC"
  } else {
    prior <- "Corr"
  }

  if (outcomeClass == "character" & !is.list(Koutter)){
    size.m <- min(table(metadata[, var2predict]))
    Kint.m <- as.integer(size.m-(size.m/Koutter))
    if (Kint.m < 2) {Kint.m <- 2}
    if (size.m < 3) {
      stop("The smallest group has too few samples, minimum number of samples per group is 3.")
    }
    if (size.m < 10) {
      warning("The smallest group has too few samples, the models may not work properly.")
    }
    if (Koutter < 2) {stop("Koutter must be 2 or more")}
    if (Kinner < 2) {stop("Kinner must be 2 or more")}
    if (Koutter > size.m) {
      stop(paste0("Koutter must be smaller than or equal to the smallest group. Maximum Koutter is ", size.m))
    }
    if (Kinner > Kint.m) {
      stop(paste0("Not enough samples for ", Koutter,
                  " Koutter and ", Kinner, " Kinner. For ", Koutter,
                  " Koutter, Kinner must be less or equal to ",
                  Kint.m))
    }
  }

  if(is.null(positiveClass)){
    positiveClass <- levels(factor(metadata[,var2predict]))[1]
  }

  if(length(unique(expData$group))>2 & outcomeClass == "character"){
    warning("glm, ada and gamboost models are not available for multi-class models")
    models<-models[!names(models) %in% c('glm', 'ada','gamboost')]
  }
  if(length(models)==0){
    stop("No algorithm suitable. Please, check the MethodsML function")
  }

  ## 2. Koutter
  if(!is.list(Koutter)){
    #sampleSets <- unname(vapply(seq_len(Koutter), function(x){
    #caret::createDataPartition(y=expData$group, p=splitProp, list=TRUE)},
    #list(seq_len(Koutter))))
    if(!is.null(paired)){
      isPaired<-metadata[rownames(expData),paired]
    }else{
      isPaired<-NULL
    }
    sampleSets <- .makeClassBalancedFolds(y=expData$group,kfold = Koutter,
                                          repeats = 1, varType = outcomeClass,
                                          paired = isPaired)
  } else {
    sampleSets <- Koutter
  }

  ntest<-sum(unlist(lapply(sampleSets,function(n){
    nrow(expData[-as.numeric(unlist(n)),])})))

  message("Training models...")
  pb <- txtProgressBar(min = 0, max = length(sampleSets)*length(models), style = 3) # Progress bar starts

  resultNested <- lapply(sampleSets, function(x){
    training <- expData[as.numeric(unlist(x)),]
    testing <- expData[-as.numeric(unlist(x)),]

    if(!is.null(paired)){
      isPaired<-metadata[rownames(training),paired]
    }else{
      isPaired<-NULL
    }

    folds<-.makeClassBalancedFolds(y=training$group,kfold = Kinner,
                                   repeats = repeatsCV,varType = outcomeClass,
                                   paired = isPaired)

    if(!is.null(filterFeatures)){

      if(filterFeatures=="sbf"){
        filterCtrl <- caret::sbfControl(functions = NULL,method = "cv",
                                        verbose = FALSE,returnResamp = "final",
                                        index = folds,allowParallel = TRUE)
      }
      if(filterFeatures=="rfe"){
        filterCtrl <- caret::rfeControl(functions = NULL,method = "cv",
                                        verbose = FALSE,returnResamp = "final",
                                        index = folds,allowParallel = TRUE,
                                        rerank=rerank)
      }

      tmp<-as.data.frame(training[,2:ncol(training)])
      if(outcomeClass=="character"){
        y<-factor(training$group)
        if(filterFeatures=="sbf"){
          filters <- caret::sbf(x = tmp, y = factor(training$group), sbfControl = filterCtrl)
        }
        if(filterFeatures=="rfe"){
          filters <- caret::rfe(x = tmp, y = factor(training$group,),
                                rfeControl = filterCtrl,sizes = filterSizes)
        }

      }else{
        y<-as.numeric(training$group)
        if(filterFeatures=="sbf"){
          filters <- caret::sbf(x = tmp, y = training$group, sbfControl = filterCtrl)
        }
        if(filterFeatures=="rfe"){
          filters <- caret::rfe(x = tmp, y = training$group,
                                rfeControl = filterCtrl, sizes = filterSizes)
        }
      }
      #print(filters)
      #print(filters$optVariables)
      if(length(filters$optVariables)>0){
        training<-training[,c("group",filters$optVariables)]
        testing<-testing[,c("group",filters$optVariables)]
      }
    }

    my_control <- caret::trainControl(method="cv", number=Kinner,
                                      savePredictions="final",
                                      #repeats=repeatsCV,
                                      classProbs=ifelse(outcomeClass=="character",
                                                        TRUE, FALSE),
                                      index= folds,
                                      search="random"
    )
    global_args <- list(group~. , training)
    global_args[["trControl"]] <- my_control

    modelList <- lapply(models, function(m) {
      model_args <- c(global_args, m)
      if (continue_on_fail == TRUE) {
        warn <- err <- NULL
        model <- .removeOutText(
          withCallingHandlers(
            tryCatch(do.call(caret::train, model_args),
                     error = function(e) {
                       err <- conditionMessage(e)
                       NULL
                     }), warning = function(w) {
                       warn <<- append(warn, conditionMessage(w))
                       invokeRestart("muffleWarning")
                     }
          )
        )
        if (!is.null(saveLogFile) & (!is.null(warn)|!is.null(err))) {
          cat(paste0("Model ", m$method, ":\n", paste0(err, collapse = "\n"),
                     paste0(warn, collapse = "\n"), "\n\n"), file = saveLogFile, append = TRUE)
        }
      } else {
        model <- .removeOutText(do.call(caret::train, model_args))
      }
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      return(model)
    })

    names(modelList) <- names(models)
    nulls <- sapply(modelList, is.null)
    modelList <- modelList[!nulls]
    if (length(modelList) == 0) {
      stop("caret:train failed for all models.  Please inspect your data.")
    }
    class(modelList) <- c("caretList")
    modelResults <- .removeOutText(modelList)

    ## Get model stats for Koutter
    predictionTable <- list()
    cm <- list()
    for (m in 1:length(modelResults)) {
      if (outcomeClass == "character") {
        classLabels <- levels(as.factor(training$group))
        predTest <- stats::predict(modelResults[[m]], newdata = testing,
                                   type = "prob")[, classLabels]
        sel <- unlist(lapply(1:nrow(predTest), function(n) {
          ifelse(sum(!is.na(predTest[n, ])) == 0, FALSE,
                 TRUE)
        }))
        if (!any(sel)) {
          predTest <- NULL
          cm <- append(cm, list(NULL))
          names(cm)[m] <- names(modelResults)[m]
        } else {
          y <- factor(testing$group[sel], levels = unique(testing$group))
          predTest <- predTest[sel, ]
          x <- factor(unlist(lapply(1:nrow(predTest),
                                    function(n) {
                                      names(which.max(predTest[n, ]))
                                    })), levels = unique(testing$group))
          cmModel <- caret::confusionMatrix(x, y, positive = positiveClass)
          cm <- append(cm, list(cmModel))
          names(cm)[m] <- names(modelResults)[m]
          predTest <- data.frame(predTest, obs = y)
        }
      } else {
        predTest <- stats::predict(modelResults[[m]], newdata = testing)
        predTest <- data.frame(pred = as.numeric(predTest),
                               obs = as.numeric(testing$group))
        cm <- NULL
      }
      predictionTable <- append(predictionTable, list(stats::na.omit(predTest)))
      names(predictionTable)[m] <- names(modelResults)[m]
    }
    gc()
    return(list(models=modelResults, preds=predictionTable, cm=cm))
  })
  close(pb) # Progress bar ends
  message("Done")

  Finalfeatures<-unique(c(unlist(lapply(resultNested,function(f.ns){
    tmp.ff<-colnames(f.ns$models[[1]]$trainingData)
    return(tmp.ff[!tmp.ff %in% ".outcome"])
  }))))

  for (s in 1:length(resultNested)) {
    m_toremove <- c()
    for (m in 1:length(resultNested[[s]][["models"]])) {
      if (is.null(resultNested[[s]][["models"]][[m]]) |
          is.null(resultNested[[s]][["preds"]][[m]])) {
        m_toremove <- append(m_toremove, m)
      }
    }
    if(!is.null(m_toremove)){
      resultNested[[s]][["models"]][[m_toremove]] <- NULL
      resultNested[[s]][["preds"]][[m_toremove]] <- NULL
      resultNested[[s]][["cm"]][[m_toremove]] <- NULL
    }
  }

  validModels <- list()
  validModels <- lapply(resultNested, function(m) {
    vm <- append(validModels, names(m$models))
    return(vm)
  })
  validModels <- unique(unlist(validModels))
  nullmodels <- do.call("cbind", lapply(1:length(resultNested), function(x) {
    nullmodel <- data.frame(row.names = validModels, !validModels%in%names(resultNested[[x]][["models"]]))
    nullmodel[!nullmodel[,1],] <- 0
    colnames(nullmodel) <- paste0("sampleset", x)
    return(nullmodel)
  }))
  perc_Koutter_nullModel <- 50 # percentage of Koutter sample sets with a valid model required to keep that model in results.
  removeModel <- rowSums(nullmodels)>length(resultNested)*(1-(perc_Koutter_nullModel/100))
  removeModel <- names(removeModel[removeModel])
  for (i in 1:length(resultNested)) {
    resultNested[[i]][["models"]][names(resultNested[[i]][["models"]])%in%removeModel] <- NULL
    resultNested[[i]][["preds"]][names(resultNested[[i]][["preds"]])%in%removeModel] <- NULL
  }
  validModels <- validModels[!validModels %in% removeModel]
  failedModels <- names(models)[!names(models) %in% validModels]
  if (length(failedModels) > 0) {
    message(paste0("The following models failed: ", paste0(failedModels,
                                                           collapse = ", ")))
  }
  if (is.null(validModels)) {
    stop("All models failed. Please inspect your data.")
  }
  models <- models[validModels]
  if (outcomeClass == "character") {
    metrics <- c("mcc", "balacc", "accuracy", "recall",
                 "specificity", "npv", "precision", "fscore")
    type <- "classification"
    levels <- c(positiveClass, unique(unique(expData$group))[!unique(expData$group) %in%
                                                               positiveClass])
  } else {
    metrics <- c("r", "RMSE", "R2", "MAE", "RMAE", "RSE")
    type <- "regression"
    levels <- NULL
  }
  message("Calculating metrics...")
  stats <- do.call("cbind", lapply(names(models), function(x) {
    sum.model.x <- do.call("cbind", lapply(1:length(sampleSets),
                                           function(it) {
                                             tmp <- resultNested[[it]]$preds[[x]]
                                             if (outcomeClass == "character" & !is.null(tmp)) {
                                               obs <- factor(tmp$obs, levels = levels)
                                               lab.pred <- unlist(lapply(1:nrow(tmp), function(n) {
                                                 names(which.max(tmp[n, !colnames(tmp) %in%
                                                                       "obs"]))
                                               }))
                                               lab.pred <- factor(lab.pred, levels = levels)
                                             } else {
                                               obs <- tmp$obs
                                               lab.pred <- tmp$pred
                                             }


                                             if (continue_on_fail == TRUE) {
                                               warn <- err <- NULL
                                               resultsTable <- withCallingHandlers(tryCatch(
                                                 metrica::metrics_summary(obs = obs,
                                                                          pred = lab.pred, type = type, pos_level = "1",
                                                                          metrics_list = metrics),
                                                 error = function(e) {
                                                   err <- conditionMessage(e)
                                                   NULL
                                                 }),
                                                 warning = function(w) {
                                                   warn <<- append(warn, conditionMessage(w))
                                                   invokeRestart("muffleWarning")
                                                 })
                                               if (!is.null(saveLogFile) & (!is.null(warn) | !is.null(err))) {
                                                 cat(paste0("Model ", x, ", metrics_summary:\n",
                                                            paste0(err, collapse = "\n"),
                                                            paste0(warn, collapse = "\n"), "\n\n"),
                                                     file = saveLogFile, append = TRUE)
                                               }
                                             } else {
                                               resultsTable <- suppressWarnings(metrica::metrics_summary(obs = obs,
                                                                                                         pred = lab.pred,
                                                                                                         type = type,
                                                                                                         pos_level = "1",
                                                                                                         metrics_list = metrics))
                                             }
                                             if (is.null(resultsTable)) {
                                               res <- data.frame(Score = rep(NA, length(metrics)))
                                             } else {
                                               rownames(resultsTable) <- resultsTable$Metric
                                               res <- data.frame(Score = resultsTable[metrics, "Score"])
                                             }
                                             rownames(res) <- metrics
                                             return(res)
                                           }))
    sum.model.all <- data.frame(results = apply(sum.model.x,
                                                1, function(metr) {
                                                  mean(metr, na.rm = T)
                                                }))
    return(sum.model.all)
  }))

  names(stats)<-names(models)

  if(length(models) > 1){
    switch(prior,
           "MCC"={stats <- stats[,order(as.numeric(stats["mcc",]),decreasing=TRUE)]},
           "Corr"={
             #rownames(stats) <- "Correlation"
             stats <- stats[,order(as.numeric(stats["r",]),decreasing=TRUE)]})
  }

  ## 4. Build model with all data, best algorithm and best parameters
  parameters <- do.call("rbind", lapply(resultNested, function(x){
    x$models[[colnames(stats)[1]]]$bestTune
  }))

  ## Table of prediction from subsets
  predsTable<-do.call("rbind",lapply(resultNested,function(x){
    x$preds[colnames(stats)[1]][[1]]
  }))

  lossSamples<-c(unlist(lapply(1:ncol(stats),function(n){
    ls.mod<-sum(unlist(lapply(resultNested,function(m){
      nrow(m$preds[colnames(stats)[n]][[1]])})))})))
  stats<-rbind(stats,"perc.lossSamples"=100-((lossSamples/ntest)*100))

  message("Done")

  bestTune <- data.frame(lapply(seq_len(ncol(parameters)), function(x){
    tmpValues <- data.frame(parameters[,x])
    if(!methods::is(tmpValues[,1], "numeric")){
      tmpValues <- names(table(tmpValues))[order(table(tmpValues), decreasing = T)][1]
    }else{
      tmpValues <- apply(tmpValues, 2, mean)
    }}))
  colnames(bestTune) <- paste0(".", colnames(parameters))
  rownames(bestTune) <- NULL

  if (".max_depth"%in%colnames(bestTune)) {bestTune$.max_depth <- round(bestTune$.max_depth)}

  rm(resultNested)
  gc()

  message("Training final model in all samples...")
  newData<-expData[,colnames(expData) %in% c("group",Finalfeatures)]

## Add aditional parameters for nnet
if (colnames(stats)[1] == "nnet" & outcomeClass == "character"){
  maxNW<-((length(Finalfeatures) + 1)*bestTune$.size[1]) +
        ((bestTune$.size[1] + 1)* length(unique(newData$group)))
  #additional_params <- list(MaxNWts = round(maxNW * 1.5,digits=0))

  fit.model <- withCallingHandlers(tryCatch(pathMED:::.removeOutText(caret::train(group ~ ., data = newData,
                                                                                  method = colnames(stats)[1],
                                                                                   tuneGrid = bestTune,
                                                                                   MaxNWts = round(maxNW * 10,digits=0))),
                                              error = function(e) {
                                                message(paste0("Error fitting the best model (", colnames(stats)[1], ") in all samples. NULL model returned. Try manually selecting a subset of samples and use the optimal parameters provided."))
                                                NULL
                                              }))

} else {

    fit.model <- withCallingHandlers(tryCatch(pathMED:::.removeOutText(caret::train(group ~ ., data = newData,
                                                                                  method = colnames(stats)[1],
                                                                                   tuneGrid = bestTune)),
                                              error = function(e) {
                                                message(paste0("Error fitting the best model (", colnames(stats)[1], ") in all samples. NULL model returned. Try manually selecting a subset of samples and use the optimal parameters provided."))
                                                NULL
                                              }))

}


  #if (continue_on_fail == TRUE) {
    #fit.model <- withCallingHandlers(tryCatch(pathMED:::.removeOutText(caret::train(group ~ ., data = newData,
    #                                                                                method = colnames(stats)[1],
    #                                                                                tuneGrid = bestTune,
    #                                                                                ... = additional_params)),
    #                                          error = function(e) {
    #                                             message(paste0("Error fitting the best model (", colnames(stats)[1], ") in all samples. NULL model returned. Try manually selecting a subset of samples and use the optimal parameters provided."))
    #                                            NULL
    #                                          }))
  #}
  message("Done")

  return(list(model=fit.model, stats=stats, bestTune=bestTune, subsample.preds = predsTable))

}

