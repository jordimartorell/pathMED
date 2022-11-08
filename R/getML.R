#' Build machine-learning based on gene-expression/mscores data
#'
#' Retrieve best model for a specific outcome
#'
#' @param expData expression/mscore matrix. Samples in columns and features in
#' rows
#' @param metadata dataframe with information for each sample. Samples in rows
#' and variables in columns
#' @param algorithms 'glm','lm','lda','xgbTree','rf','knn','svmLinear','nnet',
#' 'svmRadial','nb','lars','rpart' or 'all' (all algorithms are used)
#' @param add additional model list (supported by caret) in format:
#' add=list(modelName = caretModelSpec(method='modelName',
#' tuneGrid=(.parameters='values')))
#' @param var2predict character with the column name of the @metadata to predict
#' @param outerfolds number of external folds in which exp.data is divided
#' @param splitProp proportion of samples used as train for the outerfolds (0-1)
#' @param innerfolds number of internal folds (for parameter tuning)
#' @param innerRepeats number of repetitions of the parameter tuning process
#' @param featureFilter method to reduce number of features (none,fcbf)
#' @param prior rank best model based on AUC (AUC), the mean of AUC and the
#' balanced accuracy (BAUC), the balanced accuracy (BA) or manual (Manual)
#' @param cores Number of cores to be used.
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
#' data(exampleRefMScore, exampleMetadata)
#' \donttest{
#' relevantPaths <- diseasePaths(MRef=exampleRefMScore,
#' min_datasets=3,
#' perc_samples=10)
#'
#' MScoresExample <- GetMscores(genesets = relevantPaths,
#' Patient = exampleData,
#' Healthy = NULL,
#' nk = 5)
#'
#' fit.model <- getML(expData=MScoresExample,
#' metadata=exampleMetadata,
#' var2predict="Response",
#' algorithms="svmLinear",
#' outerfolds=10,
#' splitProp=0.8,
#' innerfolds=10,
#' innerRepeats=10,
#' featureFilter = "none",
#' prior="BAUC")
#' }
#'
#' @export
getML <- function(expData,
                  metadata,
                  algorithms='all',
                  add=NULL,
                  var2predict,
                  outerfolds=4,
                  splitProp=0.75,
                  innerfolds=10,
                  innerRepeats=5,
                  featureFilter='none',
                  prior='AUC',
                  cores=1){

    if (!var2predict %in% colnames(metadata)) {
        stop("var2predict must be a column of metadata")
    }

    ## 1. Formatting data input
    samples <- intersect(rownames(metadata), colnames(expData))

    if (length(samples) < 1) {
        stop("Row names of metadata and column names of expData do not match")
    }

    expData <- expData[,samples]
    metadata <- metadata[samples,]

    expData <- data.frame('group'=metadata[,var2predict],
                          as.data.frame(t(expData)))
    colnames(expData) <- stringi::stri_replace_all_regex(
        colnames(expData),
        pattern = c('/',' ','-'),
        replacement = c('.','.','.'),
        vectorize=FALSE)
    expData <- expData[!is.na(expData$group),]
    outcomeClass <- class(expData$group)
    if (!methods::is(expData$group, "character")){
        prior <- "Corr"
    }

    ## Get algorithm variables
    methodList <- .methodsML(algorithms=algorithms, add=NULL,
                             outcomeClass=outcomeClass, training=expData)


    ## 2. Outerfolds
    outerSplits <- unname(vapply(seq_len(outerfolds), function(x){
        createDataPartition(y=expData$group, p=splitProp, list=TRUE)},
        list(seq_len(outerfolds))))

    resOuter <- BiocParallel::bplapply(outerSplits, function(x){
        training <- expData[as.numeric(unlist(x)),]
        testing <- expData[-as.numeric(unlist(x)),]
        my_control <- trainControl(method="repeatedcv", number=innerfolds,
                                   savePredictions="final",
                                   repeats=innerRepeats,
                                   classProbs=ifelse(outcomeClass=="character",
                                                     TRUE, FALSE),
                                   index=createResample(training$group,
                                                        innerfolds))

        invisible(switch(featureFilter,
                         "fcbf"={
                             training <- .fast.cor.FS(training, thresh=0.0025)
                             testing <- testing[,colnames(training)]},
                         "none"={}
        ))

        modelResults <- caretList(group~., data=training,
                                  trControl=my_control,
                                  tuneList=methodList)

        ## Get model stats for outerfolds
        predictionTable <- list()
        cm <- list()
        for(model in modelResults){
            if(outcomeClass=="character"){ ## Categorical outcome
                classLabels <- levels(as.factor(training$group))
                predTest <- stats::predict(model, newdata=testing,
                                           type="prob")[,classLabels]
                x <- as.factor(colnames(predTest)[apply(predTest, 1,
                                                        which.max)])
                y <- as.factor(testing$group)
                cmModel <- confusionMatrix(x,y)
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
    }, BPPARAM=BiocParallel::SnowParam(workers = cores, progressbar=TRUE))

    ## 3. Best algorithm selection (model prioritization)
    stats <- as.data.frame(do.call("cbind",
                                   lapply(names(methodList),
      function(mod){
          if(outcomeClass=="character"){
              cm.mod <- do.call("cbind", lapply(resOuter,function(x){
                  c(x$cm[[mod]]$overall["Accuracy"], x$cm[[mod]]$byClass)
              }))

              allpreds <- do.call("rbind", lapply(resOuter,
                                                  function(x){x$preds[[mod]]}))
              AUC <- invisible(MLeval::evalm(allpreds))$stdres[[1]]["AUC-ROC",
                                                                    "Score"]
              statsTmp <- c(AUC, apply(cm.mod, 1, mean))
              names(statsTmp)[1] <- "AUC"
              return(statsTmp)

          } else{
              allpreds <- do.call("rbind", lapply(resOuter,
                                                  function(x){x$preds[[mod]]}))
              statsTmp <- stats::cor(allpreds$pred, allpreds$obs)

          }
      })))

    colnames(stats) <- names(methodList)

    if(length(methodList) > 1){
        switch(prior,
               "AUC"={stats <- stats[,order(stats["AUC",],decreasing=TRUE)]},
               "BAUC"={
                   ord <- apply(stats[c("AUC","Balanced Accuracy"),], 2 ,sum)
                   stats <- stats[,order(ord,decreasing=TRUE)]},
               "Manual"={
                   print(stats)
                   sel <- as.character(readline("Select algorithm > "))
                   stats <- stats[,c(sel,setdiff(colnames(stats),sel))]},
               "BA"={stats <- stats[,order(as.numeric(
                   stats["Balanced Accuracy",]),
                   decreasing=TRUE)]},
               "Corr"={
                   rownames(stats) <- "Correlation"
                   stats <- stats[,order(as.numeric(stats["Correlation",]),
                                         decreasing=TRUE)]})
    }

    ## 4. Build model with all data, best algorithm and best parameters
    parameters <- do.call("rbind", lapply(resOuter, function(x){
        x$models[[colnames(stats)[1]]]$bestTune
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

    invisible(switch(featureFilter,
                     "fcbf"={
                         expData <- .fast.cor.FS(expData, thresh=0.0025)},
                     "none"={}
    ))

    my_control <- trainControl(method="repeatedcv", number=innerfolds,
                               savePredictions="final", repeats=innerRepeats,
                               classProbs=ifelse(outcomeClass=="character",
                                                 TRUE, FALSE))

    fit.model <- train(group~.,data=expData, method=colnames(stats)[1],
                       tuneGrid=bestTune, trControl=my_control)
    auc <- NULL
    if(outcomeClass=="character"){
        allpreds <- do.call("rbind",
                            lapply(resOuter,
                                   function(x){x$preds[[colnames(stats)[1]]]}))
        auc <- MLeval::evalm(allpreds)
    }
    return(list(model=fit.model, stats=stats, bestTune=bestTune, auc=auc))
}


