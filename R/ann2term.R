#' Annotate pathways
#'
#' @param MScores Matrix with pathways IDs as row names
#'
#' @return A dataframe with the input IDs and their corresponding annotations
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{getMscoresRef}}
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
#' }
#' MScoresAnnotated <- ann2term(MScoresExample)
#' @export
ann2term <- function(MScores){
    '%>%' <- magrittr::`%>%`
    MScores_df <- as.data.frame(MScores)
    anns <- rownames(MScores)
    splits <- sapply(strsplit(anns,split = ".split"),"[",2)
    splits <- unlist(lapply(splits, function(x){
        if (!is.na(x)){
            x = paste0(".split",x)
        } else{
            x = ""
        }
        return(x)
    }))
    anns_prev <- sapply(strsplit(anns,split = ".split"),"[",1)
    anns_df <- data.frame(ann_prev = anns_prev, split = splits, full = paste0(anns_prev,splits))
    anns_merge <- ann_info %>% dplyr::inner_join(anns_df, by = c("annotation_id" = "ann_prev")) %>% as.data.frame()
    rownames(anns_merge) <- anns_merge$full
    anns_merge <- anns_merge[anns,]
    anns_merge <- paste0(anns_merge$term,anns_merge$split)
    tableOut <- data.frame(ID = rownames(MScores), term = anns_merge)
    return(tableOut)
}
