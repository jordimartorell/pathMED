#' Prepare the models parameter for the trainModel function
#'
#' @importFrom caretEnsemble caretModelSpec
#' @param algorithms Vector with one or more of these methods: 'glm', 'lm',
#' 'lda', 'xgbTree', 'rf', 'knn', 'svmLinear', 'nnet', 'svmRadial', 'nb',
#' 'lars','rpart', 'gamboost', 'ada', 'brnn', 'enet', or 'all' to use all
#' algorithms
#' @param outcomeClass Predicted variable type ('character' or 'numeric')
#' @param tuneLength maximum number of tuning parameter combinations
#'
#' @return A list with the selected models ready to use as the 'models'
#' parameter in the trainModel function
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{trainModel}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' models <- methodsML(c("rf", "knn"), tuneLength = 20,
#'                     outcomeClass = "character")
#' @export
methodsML <- function(
        algorithms = c("rf", "knn", "nb"),
        outcomeClass,
        tuneLength = 20) {
    if (is.null(outcomeClass)) {
        stop("outcomeClass is missing, please specify whether the variable to
            predict is character or numeric.")
    }
    if ("all" %in% algorithms) {
        algorithms <- c(
            "glm", "lm", "lda", "xgbTree", "rf", "knn", "svmLinear",
            "svmRadial", "nnet", "nb", "lars", "rpart", "ada",
            "gamboost", "brnn", "enet"
        )
    }
    availableMethods <- list(
        numeric = c(
            "glm", "lm", "xgbTree", "rf", "knn",
            "nnet", "svmRadial", "svmLinear", "lars",
            "rpart", "gamboost", "brnn", "enet"
        ),
        character = c(
            "glm", "lda", "xgbTree", "rf", "knn",
            "nnet", "svmLinear", "svmRadial", "nb",
            "ada", "gamboost"
        )
    )

    metric <- ifelse(outcomeClass == "character", "Kappa", "RMSE")

    methodList <- list(
        lm = caretModelSpec(
            method = "lm", tuneLength = tuneLength,
            metric = metric
        ),
        glm = caretModelSpec(
            method = "glm", tuneLength = tuneLength,
            metric = metric
        ),
        lda = caretModelSpec(
            method = "lda", tuneLength = tuneLength,
            metric = metric
        ),
        xgbTree = caretModelSpec(
            method = "xgbTree", tuneLength = tuneLength,
            metric = metric, nthread = 1
        ),
        rf = caretModelSpec(
            method = "rf", tuneLength = tuneLength,
            metric = metric
        ),
        knn = caretModelSpec(
            method = "knn", tuneLength = tuneLength,
            metric = metric
        ),
        svmLinear = caretModelSpec(
            method = "svmLinear", tuneLength = tuneLength,
            metric = metric
        ),
        svmRadial = caretModelSpec(
            method = "svmRadial", tuneLength = tuneLength,
            metric = metric
        ),
        nnet = caretModelSpec(
            method = "nnet", tuneLength = tuneLength,
            metric = metric
        ),
        nb = caretModelSpec(
            method = "nb", tuneLength = tuneLength,
            metric = metric
        ),
        lars = caretModelSpec(
            method = "lars", tuneLength = tuneLength,
            metric = metric
        ),
        rpart = caretModelSpec(
            method = "rpart", tuneLength = tuneLength,
            metric = metric
        ),
        ada = caretModelSpec(
            method = "ada", tuneLength = tuneLength,
            metric = metric
        ),
        gamboost = caretModelSpec(
            method = "gamboost", tuneLength = tuneLength,
            metric = metric
        ),
        brnn = caretModelSpec(
            method = "brnn", tuneLength = tuneLength,
            metric = metric
        ),
        enet = caretModelSpec(
            method = "enet", tuneLength = tuneLength,
            metric = metric
        )
    )


    algorithmsGood <- algorithms[algorithms %in% unlist(
        availableMethods[outcomeClass]
    )]

    if (length(algorithmsGood) < 1) {
        stop(
            "Selected algorithms do not work with", outcomeClass,
            "outcome variables."
        )
    }

    removedAlg <- setdiff(algorithms, algorithmsGood)

    if (length(removedAlg) > 0 & !("all" %in% algorithms)) {
        warning(
            "Methods", paste(removedAlg, collapse = ", "),
            "have been removed. Not suitable for", outcomeClass,
            "outcome variables."
        )
    }

    methodList <- methodList[algorithmsGood]

    return(methodList)
}
