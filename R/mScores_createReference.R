#' Create a reference dataset based on M-scores
#'
#' @param refObject A refObject object structure: a list of lists, each one
#' with a
#' cases expression matrix and controls expression matrix
#' (named as Disease and Healthy). It can be constructed with the buildRefObject
#'  function.
#' @param geneSets A named list with each gene set,
#' or the name of one preloaded database (go_bp, go_cc, go_mf, kegg, reactome,
#' pharmgkb, lincs, ctd, disgenet, hpo, wikipathways, tmod)
#' or a GeneSetCollection.
#' @param cores Number of cores to be used.
#'
#' @return A list with three elements. The first one is a list with the M-scores
#' for each dataset. The second one is the geneSet used for the analysis and
#' the third one is the input data.
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{mScores_imputeFromReference}}, \code{\link{dissectDB}},
#' \code{\link{mScores_filterPaths}}, \code{\link{trainModel}}
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
#' refObject <- buildRefObject(
#'     data = list(
#'         refData$dataset1, refData$dataset2,
#'         refData$dataset3, refData$dataset4
#'     ),
#'     metadata = list(
#'         refData$metadata1, refData$metadata2,
#'         refData$metadata3, refData$metadata4
#'     ),
#'     groupVar = "group",
#'     controlGroup = "Healthy_sample"
#' )
#'
#' refMscore <- mScores_createReference(refObject, geneSets = "tmod")
#' @export
mScores_createReference <- function(
        refObject,
        geneSets,
        cores = 1) {

    if (is(geneSets, "GeneSetCollection")) {
        geneSets <- .gsc_to_list(geneSets)
    }
    else if (!is(geneSets, "list")) {
        geneSets <- genesetsData[[geneSets]]
    }

    nDatasets <- length(refObject)

    mscores <- lapply(seq_len(nDatasets), function(i) {
        dataset <- refObject[[i]]
        PatientData <- as.data.frame(dataset[[1]])
        HealthyData <- as.data.frame(dataset[[2]])
        HealthyData <- HealthyData[apply(HealthyData, 1, stats::sd) != 0, ]
        HealthyData <- HealthyData[!is.na(rownames(HealthyData)), ]
        HealthyMean <- rowMeans(HealthyData, na.rm = TRUE)
        HealthySD <- apply(HealthyData, 1, function(x) {
            sd(x,
                na.rm = TRUE
            )
        })
        HealthyMeanSD <- cbind(HealthyMean, HealthySD)
        HealthyMeanSD <- HealthyMeanSD[HealthyMeanSD[, 2] != 0, ]
        message("Running dataset ", i, " of ", nDatasets)

        # Avoid using more cores than samples
        workers <- min(cores, ncol(PatientData))

        res <- BiocParallel::bplapply(seq_len(ncol(PatientData)),
            function(
        column, PatientData, geneNames,
        geneSets, HealthyMeanSD,
        .getMscorePath) {
                Patient <- PatientData[, column]
                names(Patient) <- geneNames
                res.i <- lapply(geneSets,
                    .getMscorePath,
                    HealthyMeanSD = HealthyMeanSD,
                    Patient = Patient
                )
                res.i <- as.data.frame(do.call("rbind", res.i))
                return(res.i)
            },
            PatientData = PatientData,
            geneNames = rownames(PatientData),
            geneSets = geneSets,
            HealthyMeanSD = HealthyMeanSD,
            .getMscorePath = .getMscorePath,
            BPPARAM = BiocParallel::SnowParam(
                workers = workers, progressbar = TRUE,
                exportglobals = FALSE
            )
        )

        res <- do.call("cbind", res)
        colnames(res) <- colnames(PatientData)
        return(res)
    })
    if (!is.null(names(refObject))) {
        names(mscores) <- names(refObject)
    }
    return(list(mscores = mscores, geneSets = geneSets, input = refObject))
}
