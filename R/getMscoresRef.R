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
#' data(refData)
#' refMscore <- getMscoresRef(data=refData, genesets="tmod")
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

    lengthData <- length(data)

    data.Mscore <- lapply(seq_len(lengthData), function (i) {
        dataset <- data[[i]]
        Patient <- dataset[[1]]
        Healthy <- as.matrix(dataset[[2]])
        Healthy <- Healthy[ifelse(apply(Healthy, 1, stats::sd) == 0, FALSE,
                                  TRUE),]
        Healthy <- Healthy[!is.na(rownames(Healthy)),]
        H <- data.frame(rowMeans(Healthy),
                      matrixStats::rowSds(Healthy))
        rownames(H) <- rownames(Healthy)
        message("Running dataset ", i, " of ", lengthData)
        res <- BiocParallel::bplapply(Patient, function(pat, geneNames,
                                                        path.list, Healthy) {
            names(pat) <- geneNames
            res.i <- lapply(path.list,
                            .getMscorePath,
                            Healthy=Healthy,
                            Patient=pat)
            res.i <- as.data.frame(do.call("rbind", res.i))
            return(res.i)
        },
        geneNames=rownames(Patient),
        path.list=path.list,
        Healthy=H,
        BPPARAM=BiocParallel::SnowParam(workers=cores, progressbar=TRUE))

        res <- do.call("cbind", res)
        colnames(res) <- colnames(Patient)
        return(res)
    })
    return(list(mscores=data.Mscore, genesets=path.list, expression=data))
}
