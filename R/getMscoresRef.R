#' Calculate M-scores from reference data
#'
#' Retrieve Mscore for a list of datasets
#'
#' @param data A list of lists, each one with a cases expression matrix and
#' controls expression matrix, always in this order.
#' @param genesets character, name of the preloaded database (gobp, gomf, gocc,
#' kegg, reactome, tmod) or 'custom' to provide the annotation.
#' @param customGeneset Only if genesets == 'custom'. A named list with each
#' gene set.
#' @param cores Number of cores to be used.
#'
#' @return A list with three elements. The first one is a list with the M-scores
#' for each dataset. The second one is the geneset used for the analysis and
#' the third one is the input gene expression data.
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{diseasePaths}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' library(pathMED)
#' data(refData)
#' \donttest{
#' refMscore <- getMscoresRef(data=refData, genesets="tmod")
#' }
#' @export
getMscoresRef <- function(data,
                          genesets = "reactome",
                          customGeneset=NULL,
                          cores = 1){
    ## Create path.list
    if (genesets == "custom"){
        path.list <- customGeneset
    } else {
        path.list <- genesetsData[[genesets]]
    }

    data.Mscore <- pbapply::pblapply(data, function (dataset) {
        Patient <- dataset[[1]]
        Healthy <- dataset[[2]]

        res <- apply(Patient, 2, function(pat) {
            names(pat) <- rownames(Patient)
            res.i <- BiocParallel::bplapply(path.list,
                                            .getMscorePath,
                                            Healthy=Healthy,
                                            Patient=pat,
                                            BPPARAM=
                                                BiocParallel::SnowParam(
                                                workers = cores))
            res.i <- as.data.frame(do.call("rbind", res.i))
            return(res.i)
        })
        res <- do.call("cbind", res)
        colnames(res) <- colnames(Patient)
        return(res)
    })

    return(list(mscores=data.Mscore, genesets=path.list, expression=data))
}
