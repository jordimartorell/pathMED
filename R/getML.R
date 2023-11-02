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
#' @param Koutter number of Outter folds in which exp.data is split for neasted cold validation.
#' A list of integer with elements for each resampling iteration is admitted. 
#' Each list element is a vector of integers corresponding to the rows used for training at that iteration.
#' @param Kinner number of cross-validation folds (for parameter tuning)
#' @param repeatsCV number of repetitions of the parameter tuning process
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
#' exampleRefMScore <- getMscoresRef(data=refData, genesets="tmod")
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
                  models=methodsML(),
                  var2predict,
                  Koutter=5,
                  Kinner=4,
                  repeatsCV=5,
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

    ## 1. Formatting data input
    samples <- intersect(rownames(metadata), colnames(expData))

    if (length(samples) < 1) {
        stop("Row names of metadata and column names of expData do not match")
    }

    expData <- expData[,samples]
    metadata <- metadata[samples,,drop=FALSE]

    expData <- data.frame('group'=metadata[,var2predict],
                          as.data.frame(t(expData)))
    colnames(expData) <- stringi::stri_replace_all_regex(
        colnames(expData),
        pattern = c('/',' ','-'),
        replacement = c('.','.','.'),
        vectorize=FALSE)
    expData <- expData[!is.na(expData$group),]
    if(is.factor(expData$group)){
        expData$group<-as.character(expData$group)
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

    if(length(unique(expData$group))>2){
      warning("glm, ada and gamboost models are not available for multi-class models")
      models<-models[!names(models) %in% c('glm', 'ada','gamboost')]
      if(length(models)==0){
        stop("Any algorithm suitable. Please, check the MethodsML function")
      }
    }
  
    ## 2. Koutter
    if(!is.list(Koutter)){
      #sampleSets <- unname(vapply(seq_len(Koutter), function(x){
        #caret::createDataPartition(y=expData$group, p=splitProp, list=TRUE)},
        #list(seq_len(Koutter))))
      sampleSets <- .makeClassBalancedFolds(y=expData$group,kfold = Koutter,
                              repeats = 1, varType = outcomeClass)
        
    } else {
    sampleSets <- Koutter
    }
  
  message("Training models...")
  pb <- txtProgressBar(min = 0, max = length(sampleSets)*length(models), style = 3) # Progress bar starts
  
    resultNested <- lapply(sampleSets, function(x){
        training <- expData[as.numeric(unlist(x)),]
        testing <- expData[-as.numeric(unlist(x)),]
      
        folds<-.makeClassBalancedFolds(y=training$group,kfold = Kinner,
                              repeats = repeatsCV,varType = outcomeClass)
      
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
            model <- tryCatch(do.call(train, model_args), 
                          error = function(e) {
                            if (!is.null(saveLogFile)) {
                              cat(paste0("Error in model ", m$method, ":\n", paste0(e$message, collapse = "\n"), "\n\n"), file = saveLogFile, append = TRUE)
                            }
                            e <- NULL
                          }, 
                          warning = function(w) {
                            if (!is.null(saveLogFile)) {
                              cat(paste0("Warning in model ", m$method, ":\n", paste0(w$message, collapse = "\n"), "\n\n"), file = saveLogFile, append = TRUE)
                            }
                            w <- NULL
                          }
            )
          } else {
            model <- do.call(train, model_args)
          }
          setTxtProgressBar(pb, getTxtProgressBar(pb) + 1) # Progress bar increases
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
        for(model in modelResults){
            if(outcomeClass=="character"){ ## Categorical outcome
                classLabels <- levels(as.factor(training$group))
                predTest <- stats::predict(model, newdata=testing,
                                           type="prob")[,classLabels]

                 sel<-unlist(lapply(1:nrow(predTest),function(n){
                  ifelse(sum(!is.na(predTest[n,]))==0,FALSE,TRUE)}))
                 y <- as.factor(testing$group[sel])
                 predTest<-predTest[sel,]
                 x<-as.factor(unlist(lapply(1:nrow(predTest),function(n){
                  names(which.max(predTest[n,]))})))
              
                cmModel <- confusionMatrix(x,y,positive = positiveClass)
                cm <- append(cm,list(cmModel))
                predTest <- data.frame(predTest,"obs"=y)

            } else{ ## Numeric outcome
                predTest <- stats::predict(model, newdata=testing)
                predTest <- data.frame("pred"=as.numeric(predTest),
                                       "obs"=as.numeric(testing$group))
            }
            predictionTable <- append(predictionTable,
                                      list(stats::na.omit(predTest)))
        }

        names(predictionTable) <- names(modelResults)
        if(!!length(cm)){
            names(cm) <- names(modelResults)
        }else{
            cm <- NULL
        }
        return(list(models=modelResults, preds=predictionTable, cm=cm))
    })
          
    close(pb) # Progress bar ends
    message("Done")

    validModels <- list()
    validModels <- lapply(resultNested, function(m){
      vm <- append(validModels, names(m$models))
      return(vm)
      })
    validModels <- unique(unlist(validModels))
    failedModels <- names(models)[!names(models)%in%validModels]
    message(paste0("The following models failed: ", paste0(failedModels, collapse = ", ")))
    if (!is.null(saveLogFile)) {
      message(paste0("Error and warning messages saved in ", saveLogFile))
    }
  
    models <- models[validModels] # Remove failed models
    
    ## 3. Best algorithm selection (model prioritization)
   if(outcomeClass=="character"){
      metrics<-c("mcc","balacc","accuracy","recall","specificity","npv",
               "precision","fscore")
      type<-"classification"
      levels<-c(positiveClass,unique(unique(expData$group))[!unique(expData$group) %in% positiveClass])
  }else{
      metrics<-c("r","RMSE","R2","MAE","RMAE","RSE")
      type<-"regression"
      levels <- NULL
  }
  
  
  stats<-do.call("cbind",lapply(names(models),function(x){
    
    sum.model.x<-do.call("cbind",lapply(1:length(sampleSets),function(it){
      
      tmp<-resultNested[[it]]$preds[[x]]
      if(outcomeClass=="character"){
        obs<-factor(tmp$obs,levels =levels)
        lab.pred<-unlist(lapply(1:nrow(tmp),function(n){
          names(which.max(tmp[n,!colnames(tmp) %in% "obs"]))}))
        lab.pred<-factor(lab.pred,levels = levels)
      }else{
        obs<-tmp$obs
        lab.pred<-tmp$pred
      }
      
      resultsTable <- suppressWarnings(metrica::metrics_summary(obs = obs,
                                                                pred = lab.pred, type=type, 
                                                                pos_level = "1",
                                                                metrics_list = metrics))
      rownames(resultsTable)<-resultsTable$Metric
      res<-data.frame("Score"=resultsTable[metrics,"Score"]) 
      rownames(res)<-metrics
      return(res)
    }))
    
    sum.model.all<-data.frame("results"=apply(sum.model.x,1,function(metr){
      mean(metr,na.rm=T)}))
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

    bestTune <- data.frame(lapply(seq_len(ncol(parameters)), function(x){
        tmpValues <- data.frame(parameters[,x])
        if(!methods::is(tmpValues[,1], "numeric")){
            tmpValues <- names(table(tmpValues))[order(table(tmpValues))][1]
        }else{
            tmpValues <- apply(tmpValues, 2, mean)
        }}))
    colnames(bestTune) <- paste0(".", colnames(parameters))
    rownames(bestTune) <- NULL
  
    if (".max_depth"%in%colnames(bestTune)) {bestTune$.max_depth <- round(bestTune$.max_depth)}

    fit.model <- .removeOutText(train(group~.,data=expData,
                                      method=colnames(stats)[1],
                                      tuneGrid=bestTune))

    return(list(model=fit.model, stats=stats, bestTune=bestTune, subsample.preds = predsTable))
}

