
#' Apply predictive model to a new dataset
#'
#' @param testData Must be a numerical matrix containing the same features used
#' for the model construction in columns, and the samples (new observations) in
#' the rows
#' @param model getML output or a caret-like model object
#' @param orig.values Optional, named vector (for numerical variables) or
#' named factor (for categorical variables) with real values for each sample
#' @param positiveClass Optional, positive class to get confusion matrix.
#' Only needed when orig.values = TRUE and for categorical variables

#' @return A dataframe with prediction for the new observations (if orig.values
#' is not provided) or a list with the dataframe with predictions and a
#' dataframe with the performance metrics (if orig.values is provided)
#'
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @import stringi
#' @import stats
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
#' @export
predictExternal<- function(testData,
                  model,
                  orig.values=NULL,
                  positiveClass=NULL){

  ## CHeck if model is contained in a getML output object of not
  if ("model" %in% names(model)){
    model<-model$model
  }

  ## Check thepresence of all features needed
  colnames(testData) <- stringi::stri_replace_all_regex(
    colnames(testData), pattern = c('/',' ','-',':'),
    replacement = c('.','.','.','.'), vectorize=FALSE)

  features<-colnames(model$trainingData)[!grepl("outcome",colnames(model$trainingData))]
  if (!all(features %in% colnames(testData))) {
    stop("Missing features in testData")
  }
  testData<-testData[,features]

  ## Get predictions
  test.predictions <- stats::predict(model, newdata = testData)
  names(test.predictions)<-rownames(testData)


  ## Return results (no orig.values provided)
  if(is.null(orig.values)){
    test.predictions<-data.frame("Obs"=names(test.predictions),
                                 "Prediction"=test.predictions)
    rownames(test.predictions)<-NULL
    return(test.predictions)

  }else{
    ## Get performance woth respect to origiganl values

    if(!is.numeric(model$trainingData[,grepl("outcome",colnames(model$trainingData))])){
      if(is.null(positiveClass)){
        positiveClass<-as.character(orig.values[1])
        cat(paste0("Positive class not provided, selected: '",positiveClass,"'\n"))
      }
      type <- "classification"
      metrics <-c("mcc", "balacc", "accuracy", "recall","specificity", "npv",
                  "precision", "fscore")

      levels = c(positiveClass, setdiff(unique(orig.values), positiveClass))
      obs <- factor(orig.values, levels = levels)
      preds<-factor(test.predictions,levels = levels)

    }else{
      type <- "regression"
      metrics <-c("r", "RMSE", "R2", "MAE", "RMAE", "RSE")

      levels = c(positiveClass, setdiff(unique(orig.values), positiveClass))
      obs <- obs
      preds<-test.predictions
    }

    stats<-metrica::metrics_summary(obs = obs, pred = preds,
                                    type = type, pos_level = "1",
                                    metrics_list = metrics)

    test.predictions<-data.frame("Obs"=names(test.predictions),
                                 "Prediction"=test.predictions)
    rownames(test.predictions)<-NULL

    res<-list("Predictions"=test.predictions,
              "stats"=stats)
    return(res)
  }
}
