#' Annotate the pathways from a scores matrix
#'
#' @param scoresMatrix Matrix with pathways IDs as row names
#'
#' @return A data frame with the input IDs and their corresponding terms
#'
#' @author Raúl López-Domínguez, \email{raul.lopez@@genyo.es}
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{getScores}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(exampleData)
#' scoresExample <- getScores(exampleData, geneSets = "tmod", method = "GSVA")
#' annotatedTerms <- ann2term(scoresExample)
#' @export
ann2term <- function(scoresMatrix) {
    scoresMatrix_df <- as.data.frame(scoresMatrix)
    anns <- rownames(scoresMatrix)
    splits <- vapply(strsplit(anns, split = ".split"), "[", 2,
        FUN.VALUE = character(1)
    )
    splits <- unlist(lapply(splits, function(x) {
        if (!is.na(x)) {
            x <- paste0(".split", x)
        } else {
            x <- ""
        }
        return(x)
    }))
    anns_prev <- vapply(strsplit(anns, split = ".split"), "[", 1,
        FUN.VALUE = character(1)
    )
    anns_df <- data.frame(
        ann_prev = anns_prev, split = splits,
        full = paste0(anns_prev, splits)
    )
    anns_merge <- ann_info %>%
        dplyr::inner_join(anns_df,
            by = c(
                "annotation_id" =
                    "ann_prev"
            )
        ) %>%
        as.data.frame()
    rownames(anns_merge) <- anns_merge$full
    anns_merge <- anns_merge[anns, ]
    anns_merge <- paste0(anns_merge$term, anns_merge$split)
    tableOut <- data.frame(ID = rownames(scoresMatrix), term = anns_merge)
    return(tableOut)
}
