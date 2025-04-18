#' Predict conditions in external datasets
#'
#' @param testData Numerical matrix or data frame with the same features used
#' for the model construction in rows, and the samples (new observations) in
#' columns. An ExpressionSet may or SummarizedExperiment may also be used.
#' @param model trainModel output or a caret-like model object
#' @param realValues Optional, named vector (for numerical variables) or
#' named factor (for categorical variables) with real values for each sample
#' @param positiveClass Optional, positive class to get confusion matrix.
#' Only needed when realValues = TRUE and for categorical variables
#' @param use.assay If SummarizedExperiments are used, the number of the assay 
#' to extract the data.
#' 
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
#' data(refData)
#'
#' commonGenes <- intersect(rownames(refData$dataset1),
#'                          rownames(refData$dataset2))
#' dataset1 <- refData$dataset1[commonGenes, ]
#' dataset2 <- refData$dataset2[commonGenes, ]
#'
#' scoresExample <- getScores(dataset1, geneSets = "tmod", method = "Z-score")
#'
#' set.seed(123)
#' trainedModel <- trainModel(
#'     inputData = scoresExample,
#'     metadata = refData$metadata1,
#'     var2predict = "group",
#'     models = methodsML("svmLinear",
#'         outcomeClass = "character"
#'     ),
#'     Koutter = 2,
#'     Kinner = 2,
#'     repeatsCV = 1
#' )
#'
#' externalScores <- getScores(dataset2, geneSets = "tmod", method = "Z-score")
#' realValues <- refData$metadata2$group
#' names(realValues) <- rownames(refData$metadata2)
#' predictions <- predictExternal(externalScores, trainedModel,
#'     realValues = realValues
#' )
#'
#' print(predictions)
#'
#' @export
predictExternal <- function(
        testData,
        model,
        realValues = NULL,
        positiveClass = NULL,
        use.assay = 1) {

    if (is(testData, "data.frame")) {
        testData <- as.matrix(testData)
    }
    if (is(testData, "ExpressionSet")) {
        testData <- Biobase::exprs(testData)
    }
    if (is(testData, "SummarizedExperiment")) {
        testData <- as.matrix(SummarizedExperiment::assay(testData, use.assay))
    }
    ## CHeck if model is contained in a trainModel output object of not
    if ("model" %in% names(model)) {
        model <- model$model
    }

    testData <- t(testData)

    ## Check thepresence of all features needed
    colnames(testData) <- stringi::stri_replace_all_regex(
        colnames(testData),
        pattern = c("/", " ", "-", ":"),
        replacement = c(".", ".", ".", "."), vectorize = FALSE
    )

    features <- colnames(model$trainingData)[!grepl(
        "outcome",
        colnames(
            model$trainingData
        )
    )]
    if (!all(features %in% colnames(testData))) {
        missingFeatures <- setdiff(features, colnames(testData))
        stop("Missing features in testData: ", paste(missingFeatures,
            collapse = ", "
        ))
    }

    testData <- testData[, features]

    ## Get predictions
    test.predictions <- stats::predict(model, newdata = testData)
    names(test.predictions) <- rownames(testData)


    ## Return results (no real Values provided)
    if (is.null(realValues)) {
        test.predictions <- data.frame(
            "Obs" = names(test.predictions),
            "Prediction" = test.predictions
        )
        rownames(test.predictions) <- NULL
        return(test.predictions)
    } else {
        ## Get model performance

        if (!is.numeric(model$trainingData[, grepl(
            "outcome",
            colnames(
                model$trainingData
            )
        )])) {
            if (is.null(positiveClass)) {
                positiveClass <- sort(unique(realValues),
                    decreasing = TRUE
                )[1]
                message(
                    "Positive class not provided, selected: '",
                    positiveClass, "'\n"
                )
            }
            type <- "classification"
            metrics <- c(
                "mcc", "balacc", "accuracy", "recall", "specificity",
                "npv", "precision", "fscore"
            )

            levels <- c(positiveClass, setdiff(
                unique(realValues),
                positiveClass
            ))
            obs <- factor(realValues, levels = levels)
            preds <- factor(test.predictions, levels = levels)
        } else {
            type <- "regression"
            metrics <- c("r", "RMSE", "R2", "MAE", "RMAE", "RSE")

            levels <- c(positiveClass, setdiff(
                unique(realValues),
                positiveClass
            ))
            obs <- obs
            preds <- test.predictions
        }

        stats <- metrica::metrics_summary(
            obs = obs, pred = preds,
            type = type,
            pos_level = 1,
            metrics_list = metrics
        )
        if (sum(is.na(stats)) > 0) {
            message("Some metrics could not be calculated and are returned as 0.0")
            stats <- apply(stats, 2, function(x) {
                if (sum(is.na(x)) > 0) {
                    corrCol <- x
                    corrCol[is.na(corrCol)] <- 0.0
                    return(corrCol)
                } else {
                    return(x)
                }
            })
        }

        test.predictions <- data.frame(
            "Obs" = names(test.predictions),
            "Prediction" = test.predictions
        )
        rownames(test.predictions) <- NULL

        res <- list(
            "Predictions" = test.predictions,
            "stats" = stats
        )
        return(res)
    }
}
