#' Predict conditions in external datasets
#'
#' @param testData Numerical matrix with the same features used
#' for the model construction in rows, and the samples (new observations) in
#' columns
#' @param model trainModel output or a caret-like model object
#' @param realValues Optional, named vector (for numerical variables) or
#' named factor (for categorical variables) with real values for each sample
#' @param positiveClass Optional, positive class to get confusion matrix.
#' Only needed when realValues = TRUE and for categorical variables

#' @return A dataframe with predictions (if realValues
#' is not provided) or a list with the dataframe with predictions and a
#' dataframe with the performance metrics (if realValues is provided)
#'
#' @author Iván Ellson, \email{ivan.ellson.l@@gmail.com }
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(exampleData, exampleMetadata)
#' \donttest{
#'
#' scoresExample <- getScores(exampleData, geneSets="tmod", method="GSVA")
#'
#' trainedModel <- trainModel(inputData=scoresExample,
#'                             metadata=exampleMetadata,
#'                             var2predict="Response",
#'                             models=methodsML("svmLinear"))
#'
#' predictions <- predictExternal(externalData, trainedModel)
#' }
#'
#' @export
predictExternal <- function(testData,
                  model,
                  realValues=NULL,
                  positiveClass=NULL){

    ## CHeck if model is contained in a trainModel output object of not
    if ("model" %in% names(model)){
        model <- model$model
    }

    testData <- t(testData)

    ## Check thepresence of all features needed
    colnames(testData) <- stringi::stri_replace_all_regex(
        colnames(testData), pattern = c('/',' ','-',':'),
        replacement = c('.','.','.','.'), vectorize=FALSE)

    features <- colnames(model$trainingData)[!grepl("outcome",
                                                    colnames(
                                                        model$trainingData))]
    if (!all(features %in% colnames(testData))) {
        stop("Missing features in testData")
    }
    testData <- testData[,features]

    ## Get predictions
    test.predictions <- stats::predict(model, newdata=testData)
    names(test.predictions) <- rownames(testData)


    ## Return results (no real Values provided)
    if(is.null(realValues)){
        test.predictions <- data.frame("Obs"=names(test.predictions),
                                       "Prediction"=test.predictions)
        rownames(test.predictions) <- NULL
        return(test.predictions)

    } else{
        ## Get model performance

        if(!is.numeric(model$trainingData[,grepl("outcome",
                                                 colnames(
                                                     model$trainingData))])){
            if(is.null(positiveClass)){
                positiveClass <- sort(unique(realValues),
                                      decreasing=TRUE)[1]
                cat(paste0("Positive class not provided, selected: '",
                           positiveClass,"'\n"))
                positiveClassOrder <- length(unique(realValues))

            } else {
                positiveClassOrder <- which(sort(unique(realValues),
                                                 decreasing=FALSE) ==
                                                positiveClass)
            }
            type <- "classification"
            metrics <-c("mcc", "balacc", "accuracy", "recall","specificity",
                        "npv", "precision", "fscore")

            levels <- c(positiveClass, setdiff(unique(realValues),
                                               positiveClass))
            obs <- factor(realValues, levels = levels)
            preds <- factor(test.predictions, levels = levels)

        } else{
            type <- "regression"
            metrics <- c("r", "RMSE", "R2", "MAE", "RMAE", "RSE")

            levels <- c(positiveClass, setdiff(unique(realValues),
                                               positiveClass))
            obs <- obs
            preds <- test.predictions
        }

        stats <- metrica::metrics_summary(obs = obs, pred = preds,
                                          type = type,
                                          pos_level = positiveClassOrder,
                                          metrics_list = metrics)

        test.predictions <- data.frame("Obs"=names(test.predictions),
                                       "Prediction"=test.predictions)
        rownames(test.predictions) <- NULL

        res <- list("Predictions"=test.predictions,
                    "stats"=stats)
        return(res)
    }
}
