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
        Patient <- as.data.frame(dataset[[1]])
        Healthy <- as.data.frame(dataset[[2]])
        Healthy <- Healthy[ifelse(apply(Healthy, 1, stats::sd) == 0, FALSE,
                                  TRUE),]
        Healthy <- Healthy[!is.na(rownames(Healthy)),]
        H <- data.frame(apply(Healthy,1,function(x){mean(x,na.rm = T)}),
                        apply(Healthy,1,function(x){sd(x,na.rm = T)}))
        rownames(H) <- rownames(Healthy)
        H<-H[ifelse(H[,1]==0 | H[,2]==0,F,T),]
        message("Running dataset ", i, " of ", lengthData)

        # Avoid using more cores than samokes
        workers <- min(cores, ncol(Patient))

        res <- BiocParallel::bplapply(seq_len(ncol(Patient)),
                                      function(column, Patient, geneNames,
                                               path.list, Healthy, .getMscorePath) {
            pat <- Patient[,column]
            names(pat) <- geneNames
            res.i <- lapply(path.list,
                            .getMscorePath,
                            Healthy=Healthy,
                            Patient=pat)
            res.i <- as.data.frame(do.call("rbind", res.i))
            return(res.i)
        },
        Patient=Patient,
        geneNames=rownames(Patient),
        path.list=path.list,
        Healthy=H,
        .getMscorePath = .getMscorePath,
        BPPARAM=BiocParallel::SnowParam(workers=workers, progressbar=TRUE))

        res <- do.call("cbind", res)
        colnames(res) <- colnames(Patient)
        return(res)
    })
    return(list(mscores=data.Mscore, genesets=path.list, expression=data))
}
