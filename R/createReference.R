#' Create a reference of M-Scores
#'
#' @param datasetsList A list of lists, each one with a cases data matrix and
#' controls data matrix, always in this order.
#' @param geneSets A named list with each
#' gene set or the name of one preloaded database (gobp, gomf, gocc,
#' kegg, reactome, tmod)
#' @param cores Number of cores to be used.
#'
#' @return A list with three elements. The first one is a list with the M-scores
#' for each dataset. The second one is the geneset used for the analysis and
#' the third one is the input data.
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
#' refMscore <- createReference(datasetsList=refData, geneSets="tmod")
#' @export
createReference <- function(datasetsList,
                          geneSets = "reactome",
                          cores = 1){

    if(!is(geneSets, "list")) {
        geneSets <- genesetsData[[geneSets]]
    }

    nDatasets <- length(datasetsList)

    mscores <- lapply(seq_len(nDatasets), function (i) {
        dataset <- datasetsList[[i]]
        PatientData <- as.data.frame(dataset[[1]])
        HealthyData <- as.data.frame(dataset[[2]])
        HealthyData <- HealthyData[apply(HealthyData, 1, stats::sd) != 0,]
        HealthyData <- HealthyData[!is.na(rownames(HealthyData)),]
        HealthyMean <- rowMeans(HealthyData, na.rm = TRUE)
        HealthySD <- apply(HealthyData, 1, function(x) {sd(x,
                                                           na.rm = TRUE)})
        HealthyMeanSD <- cbind(HealthyMean, HealthySD)
        HealthyMeanSD <- HealthyMeanSD[HealthyMeanSD[,2] != 0,]
        message("Running dataset ", i, " of ", nDatasets)

        # Avoid using more cores than samples
        workers <- min(cores, ncol(PatientData))

        res <- BiocParallel::bplapply(seq_len(ncol(PatientData)),
                                      function(column, PatientData, geneNames,
                                               geneSets, HealthyMeanSD,
                                               .getMscorePath) {
            Patient <- PatientData[,column]
            names(Patient) <- geneNames
            res.i <- lapply(geneSets,
                            .getMscorePath,
                            HealthyMeanSD=HealthyMeanSD,
                            Patient=Patient)
            res.i <- as.data.frame(do.call("rbind", res.i))
            return(res.i)
        },
        PatientData=PatientData,
        geneNames=rownames(PatientData),
        geneSets=geneSets,
        HealthyMeanSD=HealthyMeanSD,
        .getMscorePath = .getMscorePath,
        BPPARAM=BiocParallel::SnowParam(workers=workers, progressbar=TRUE,
                                        exportglobals = FALSE))

        res <- do.call("cbind", res)
        colnames(res) <- colnames(PatientData)
        return(res)
    })
    if(!is.null(names(datasetsList))){
      names(mscores)<-names(datasetsList)
    }
    return(list(mscores=mscores, geneSets=geneSets, input=datasetsList))
}
