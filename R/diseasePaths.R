#' Find relevant pathways from the reference M-scores
#'
#' Subtitle
#'
#' @param MRef output from the getMscoresRef function
#' @param min_datasets number of datasets that each pathway must meet the
#' perc_samples threshold
#' @param perc_samples: minimun percentage of samples in a dataset in which a
#' pathway must be significant
#' @param Pcutoff: P-value cutoff for significance
#'
#' @return A list with the selected pathways
#'
#' @author Daniel Toro-Domínguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
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
#' \dontrun{
#' DATA<-readRDS(file = paste0(getwd(),sep="","/data/datasets.rds"))
#' }
#' DATA.Mscore <- GetMscoresReferencet(genesets=genesets, data=DATA)
#' head(DATA.Mscore)
#' @export
diseasePaths <- function(MRef,
                         min_datasets=round(length(MRef[[1]]) * 0.34),
                         perc_samples=10,
                         Pcutoff=0.05){
    MScores <- MRef[[1]]
    genesets <- MRef[[2]]
    expr.list <- MRef[[3]]

    HighDys.perc <- lapply(MScores, function(dat) {
        apply(dat, 1, function(x) {
            values <- as.numeric(x[!is.na(x)])
            values <- (length(values[abs(values) >=
                                        abs(qnorm(Pcutoff))])/length(x))*100
            return(values)
        })
    })
    HighDys.perc <- do.call("cbind", HighDys.perc)

    selected.path <- apply(HighDys.perc, 1, function(x){
        nperc <- length(x[x > perc_samples])
        ntimes <- ifelse(nperc > min_datasets, TRUE, FALSE)
        return(ntimes)
    })
    message("Selected paths: ", sum(selected.path))

    genesets <- genesets[selected.path]

    # Build the expression reference
    reference <- .MReference(expr.list=expr.list,
                          mscore.list = MScores,
                          geneset.list=genesets)

    return(list(genesets, reference))
}
