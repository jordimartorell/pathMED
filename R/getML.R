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
                  Koutter=4,
                  Kinner=10,
                  repeatsCV=5,
                  positiveClass=NULL){

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
      size.m <- min(table(metadata[,var2predict]))
      Kint.m <- as.integer(size.m/Koutter)

      if(size.m < 5){
        warning("The smallest group has too low number of samples, models could not work correctly")
      }

      if(Koutter > size.m){
        stop(paste0("Koutter must be less or equal to the smallest group. Max Koutter ",
            size.m))
      }    
      if(Kinner > Kint.m){
        stop(paste0("Not enough samples for ",
                  Koutter, " Koutter and ",Kinner," Kinner. For ", Koutter,
                   " Koutter, Kinner must be less or equal to ",Kint.m))
      }  
    }
  
    if(is.null(positiveClass)){
        positiveClass <- levels(factor(metadata[,var2predict]))[1]
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

    resultNested <- pbapply::pblapply(sampleSets, function(x){
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
        modelResults <- .removeOutText(caretEnsemble::caretList(group ~ ., data=training,
                                                 trControl=my_control,
                                                 tuneList=models))
        ## Get model stats for Koutter
        predictionTable <- list()
        cm <- list()
        for(model in modelResults){
            if(outcomeClass=="character"){ ## Categorical outcome
                classLabels <- levels(as.factor(training$group))
                predTest <- stats::predict(model, newdata=testing,
                                           type="prob")[,classLabels]
                x <- as.factor(ifelse(is.na(predTest[,1]) | is.na(predTest[,2]),NA,
                            ifelse(predTest[,1]==predTest[,2],NA,
                                   ifelse(predTest[,1]>predTest[,2],
                                          colnames(predTest)[1],
                                          colnames(predTest)[2]))))
                y <- as.factor(testing$group)
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

    ## 3. Best algorithm selection (model prioritization)
    stats <- as.data.frame(do.call("cbind",
                                   lapply(names(models),
                                          function(mod){
          if(outcomeClass=="character"){
              cm.mod <- do.call("cbind", lapply(resultNested,function(x){
                  c(x$cm[[mod]]$overall["Accuracy"], x$cm[[mod]]$byClass)
              }))

              allpreds <- do.call("rbind", lapply(resultNested,
                                                  function(x){x$preds[[mod]]}))
              AUCMCC <- invisible(MLeval::evalm(allpreds, silent=TRUE,
                                             showplots=FALSE,
                                             optimise="MCC"))$stdres[[
                  1]][c("MCC", "AUC-ROC"), "Score"]
              statsTmp <- c(AUCMCC, apply(cm.mod, 1, mean))
              names(statsTmp)[1:2] <- c("MCC","AUC")
              return(statsTmp)
          } else{
              allpreds <- do.call("rbind", lapply(resultNested,
                                                  function(x){
                                                      x$preds[[mod]]}))
              statsTmp <- stats::cor(allpreds$pred, allpreds$obs)

          }
          }
          )))

    colnames(stats) <- names(models)

    if(length(models) > 1){
        switch(prior,
               "MCC"={stats <- stats[,order(as.numeric(stats["MCC",]),decreasing=TRUE)]},
               "Corr"={
                   rownames(stats) <- "Correlation"
                   stats <- stats[,order(as.numeric(stats["Correlation",]),
                                         decreasing=TRUE)]})
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

    my_control <- caret::trainControl(method="repeatedcv", number=Kinner,
                               savePredictions="final", repeats=repeatsCV,
                               classProbs=ifelse(outcomeClass=="character",
                                                 TRUE, FALSE))

    fit.model <- .removeOutText(train(group~.,data=expData,
                                      method=colnames(stats)[1],
                                      tuneGrid=bestTune, trControl=my_control))

    return(list(model=fit.model, stats=stats, bestTune=bestTune, subsample.preds = predsTable))
}

