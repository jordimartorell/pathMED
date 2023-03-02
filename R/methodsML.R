#' Prepare ML models
#'
#' @param algorithms 'glm','lm','lda','xgbTree','rf','knn','svmLinear','nnet',
#' 'svmRadial','nb','lars','rpart', 'gamboost', 'ada', 'brnn', 'enet' or 'all'
#' (all algorithms are used)
#' @param outcomeClass Predicted variable type ('character' or 'numeric')
#' @param tuneLength maximum number of tuning parameter combinations
#'
#' @return A list with the selected models ready to use as the 'models'
#' parameter in the getml() function
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{getML}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' methodsML(c('rf', 'knn'), tuneLength=20)
#' @export
methodsML <- function(algorithms = "all",
                       outcomeClass = "character",
                       tuneLength = 15){
    if('all' %in% algorithms){
        algorithms <- c('glm', 'lm', 'lda', 'xgbTree', 'rf', 'knn', 'svmLinear',
                        'svmRadial', 'nnet', 'nb', 'lars', 'rpart', 'ada',
                        'gamboost', 'brnn', 'enet')
    }
    availableMethods <- list(numeric=c('glm', 'lm', 'xgbTree', 'rf', 'knn',
                                       'nnet', 'svmRadial', 'svmLinear', 'lars',
                                       'rpart','gamboost', 'brnn', 'enet'),
                             character=c('glm', 'lda', 'xgbTree', 'rf', 'knn',
                                         'nnet', 'svmLinear','svmRadial', 'nb',
                                         'ada', 'gamboost'))
    methodList <- list(
        lm=caretModelSpec(method='lm', tuneLength=tuneLength),
        glm=caretModelSpec(method='glm', tuneLength=tuneLength),
        lda=caretModelSpec(method='lda', tuneLength=tuneLength),
        xgbTree=caretModelSpec(method='xgbTree', tuneLength=tuneLength),
        rf=caretModelSpec(method='rf', tuneLength=tuneLength),
        knn=caretModelSpec(method='knn', tuneLength=tuneLength),
        svmLinear=caretModelSpec(method='svmLinear', tuneLength=tuneLength),
        svmRadial=caretModelSpec(method='svmRadial', tuneLength=tuneLength),
        nnet=caretModelSpec(method='nnet', tuneLength=tuneLength),
        nb=caretModelSpec(method='nb', tuneLength=tuneLength),
        lars=caretModelSpec(method='lars', tuneLength=tuneLength),
        rpart=caretModelSpec(method='rpart', tuneLength=tuneLength),
        ada=caretModelSpec(method='ada', tuneLength=tuneLength),
        gamboost=caretModelSpec(method='gamboost', tuneLength=tuneLength),
        brnn=caretModelSpec(method='brnn', tuneLength=tuneLength),
        enet=caretModelSpec(method='enet', tuneLength=tuneLength)
    )


    algorithmsGood <- algorithms[algorithms %in% unlist(
        availableMethods[outcomeClass])]
    if (length(algorithmsGood) < 1) {
        stop(paste0('Selected algorithms do not work with ', outcomeClass,
                    ' outcome variables.'))
    }
    removedAlg <- setdiff(algorithms, algorithmsGood)
    if (length(removedAlg) > 0 & !('all' %in% algorithms)) {
        warning(paste0('Methods ', paste(removedAlg, collapse=", "),
                       ' have been removed. Not suitable for ', outcomeClass,
                       ' outcome variables.'))
    }
    methodList <- methodList[algorithmsGood]
    return(methodList)
}
