#' Calculate pathways scores for a dataset
#'
#' @param inputData Data matrix or data frame.
#' @param geneSets A named list with each
#' gene set or the name of one preloaded database (gobp, gomf, gocc,
#' kegg, reactome, tmod).
#' @param method Scoring method: M-score, GSVA, ssGSEA, singscore, Z-score,
#' Plage, AUCell, MDT, MLM, ORA, UDT, ULM, FGSEA, norm_FGSEA, WMEAN, norm_WMEAN,
#' corr_WMEAN, WSUM, norm_WSUM or corr_WSUM.
#' @param labels (Only for M-Score) Vector with the samples class labels (0 or
#' "Healthy" for control samples). Optional.
#' @param cores Number of cores to be used.
#' @param ... Additional parameters for the scoring functions.
#'
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated. In addition, this list also contains
#' the input methData, pheno and platform objects
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{mScores_filterPaths }}, \code{\link{getML}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(exampleData)
#' MScoresExample <- getScores(exampleData, geneSets = "tmod")
#' @export
getScores <- function(inputData,
                      geneSets,
                      method = "GSVA",
                      labels = NULL,
                      cores = 1,
                      ...){

    if (is(inputData, "data.frame")) {
        inputData <- as.matrix(inputData)
    }

    if (is.null(externalReference)) {
        if (is.null(geneSets)) {
            stop("If externalReference is not provided, geneSets must be specified")
        }
        if(!is(geneSets, "list")) {
            geneSets <- genesetsData[[geneSets]]
        }
    }
    else {
        geneSets <- externalReference[["geneSets"]]
    }



    if (method %in% c("GSVA", "ssGSEA", "Z-score", "Plage")) {
        if (method == "GSVA") {
            paramMatrix <- GSVA::gsvaParam(inputData, geneSets, kcdf = "Gaussian",
                                           ...)
        }
        else if (method == "ssGSEA") {
            paramMatrix <- GSVA::ssgseaParam(inputData, geneSets, ...)
        }
        else if (method == "Z-score") {
            paramMatrix <- GSVA::zscoreParam(inputData, geneSets, ...)
        }
        else {
            paramMatrix <- GSVA::plageParam(inputData, geneSets, ...)
        }
        res <- GSVA::gsva(paramMatrix,
                          BPPARAM=BiocParallel::SnowParam(workers=cores))
    }

    else if (method == "singscore") {
        rankMatrix <- singscore::rankGenes(inputData, tiesMethod = "average")
        listSing <- lapply(geneSets, function(x) {
            singscore::simpleScore(rankData = rankMatrix, upSet = x, ...)
            })
        listScores <- sapply(listSing, function(x) x$TotalScore)
        if(is(listScores, "list")) {
            res <- do.call(rbind, listScores)
        }
        else{
            res <- t(listScores)
        }
        colnames(res) <- colnames(inputData)
    }

    else if (method %in% c("AUCell", "MDT", "MLM", "ORA", "UDT", "ULM",
                           "FGSEA", "norm_FGSEA", "WMEAN", "norm_WMEAN",
                           "corr_WMEAN", "WSUM", "norm_WSUM", "corr_WSUM")) {
        # Network necesary for these methods
        net <- do.call("rbind", lapply(seq_len(length(geneSets)), function (x) {
            res <- data.frame("source"=rep(names(geneSets)[[x]],
                                           length(geneSets[[x]])),
                              "target"=as.character(geneSets[[x]]),
                              "weight"=rep(1,length(geneSets[[x]])),
                              "mor"=rep(1,length(geneSets[[x]])))
            return(res)
        }))

        if (method %in% c("AUCell", "MDT", "MLM", "ORA", "UDT", "ULM")) {

            if (method == "AUCell") {
                scoreMatrix <- decoupleR::run_aucell(mat=inputData, network=net,
                                                     nproc=cores,
                                                     ...)
            }

            else if (method == "MDT") {
                scoreMatrix <- decoupleR::run_mdt(mat=inputData, network=net,
                                                  nproc=cores, ...)
            }

            else if (method == "MLM") {
                scoreMatrix <- decoupleR::run_mlm(mat=inputData, network=net,
                                                  ...)
            }

            else if (method == "ORA") {
                scoreMatrix <- decoupleR::run_ora(mat=inputData, network=net,
                                                  ...)
            }

            else if (method == "UDT") {
                scoreMatrix <- decoupleR::run_udt(mat=inputData, network=net,
                                                  ...)
            }

            else if (method == "ULM") {
                scoreMatrix <- decoupleR::run_ulm(mat=inputData, network=net,
                                                  ...)
            }

            scoreMatrix <- as.data.frame(scoreMatrix[,c("source","condition",
                                                        "score")])
        }

        else if (method %in% c("FGSEA", "norm_FGSEA")) {
            scoreMatrix <- decoupleR::run_fgsea(mat=inputData, network=net,
                                                nproc=cores, ...)
            if (method == "fGSEA") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="fgsea",
                                c("source","condition","score")])
            }
            else {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="norm_fgsea",
                                c("source","condition","score")])
            }
        }

        else if (method %in% c("WMEAN", "norm_WMEAN", "corr_WMEAN")) {
            scoreMatrix <- decoupleR::run_wmean(mat=inputData, network=net,
                                                ...)
            if (method == "WMEAN") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="wmean",
                                c("source","condition","score")])
            }
            else if (method == "norm_WMEAN") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="norm_wmean",
                                c("source","condition","score")])
            }
            else {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="corr_wmean",
                                c("source","condition","score")])
            }
        }

        else if (method %in% c("WSUM", "norm_WSUM", "corr_WSUM")) {
            scoreMatrix <- decoupleR::run_wsum(mat=inputData, network=net,
                                               ...)
            if (method == "WSUM") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="wsum",
                                c("source","condition","score")])
            }
            else if (method == "norm_WSUM") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="norm_wsum",
                                c("source","condition","score")])
            }
            else {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[scoreMatrix$statistic=="corr_wsum",
                                c("source","condition","score")])
            }
        }

        scoreMatrix <- reshape2::dcast(scoreMatrix, source~condition,
                                       value.var = "score")
        rownames(scoreMatrix)<-scoreMatrix$source
        scoreMatrix <- scoreMatrix[,-1]
        res <- scoreMatrix
    }

    else if (method == "M-score") {

        if (is.null(labels)) {
            stop("Labels parameter must be used
                 for M-Scores method")
        }


        if (0 %in% labels) {
            HealthyData <- inputData[, labels==0]
            PatientData <- inputData[, labels!=0]
        }
        else if ("Healthy" %in% labels) {
            HealthyData <- inputData[, labels=="Healthy"]
            PatientData <- inputData[, labels!="Healthy"]
        }
        else {
            stop("Reference samples in labels must be specified with 0 or
                 'Healthy'")
        }


        message("Healthy samples supplied. Calculating M-Scores using ",
                "healthy samples as reference for ", ncol(PatientData),
                " patients")

        HealthyMean <- rowMeans(HealthyData, na.rm = TRUE)
        HealthySD <- apply(HealthyData, 1, function(x) {sd(x,
                                                           na.rm = TRUE)})
        HealthyMeanSD <- cbind(HealthyMean, HealthySD)

        HealthyMeanSD <- HealthyMeanSD[HealthyMeanSD[,2] != 0,]

        res <- BiocParallel::bplapply(seq_len(ncol(PatientData)),
                                      function(column, geneSets,
                                               HealthyMeanSD) {
            Patient <- PatientData[,column, drop=TRUE]
            res.i <- lapply(geneSets,
                            .getMscorePath,
                            HealthyMeanSD=HealthyMeanSD,
                            Patient=Patient)
            res.i <- as.data.frame(do.call("rbind", res.i))
            return(res.i)
        },
        geneSets=geneSets,
        HealthyMeanSD=HealthyMeanSD,
        BPPARAM=BiocParallel::SnowParam(workers = cores, progressbar=TRUE))

        res <- do.call("cbind", res)
        colnames(res) <- colnames(PatientData)

    }

    return(res)
}

