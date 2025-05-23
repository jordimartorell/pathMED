#' Calculate pathways scores for a dataset
#'
#' @param inputData Matrix, data frame, ExpressionSet or SummarizedExperiment
#' with omics data. Feature names must match the gene sets nomenclature. To use
#' preloaded databases, they must be gene symbols.
#' @param geneSets A named list with each gene set,
#' or the name of one preloaded database (go_bp, go_cc, go_mf, kegg, reactome,
#' pharmgkb, lincs, ctd, disgenet, hpo, wikipathways, tmod)
#' or a GeneSetCollection. For using network methods,
#' a data frame including columns:
#' "source","target","weight" and "mor" (optional).
#' @param method Scoring method: M-Scores, GSVA, ssGSEA, singscore, Plage,
#' Z-score, AUCell, MDT, MLM, ORA, UDT, ULM, FGSEA, norm_FGSEA, WMEAN,
#' norm_WMEAN, corr_WMEAN, WSUM, norm_WSUM or corr_WSUM.
#' @param labels (Only for M-Scores) Vector with the samples class labels (0 or
#' "Healthy" for control samples). Optional.
#' @param cores Number of cores to be used.
#' @param use.assay If SummarizedExperiments are used, the number of the assay 
#' to extract the data.
#' @param ... Additional parameters for the scoring functions.
#'
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated. In addition, this list also contains
#' the input methData, pheno and platform objects
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{trainModel}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(pathMEDExampleData)
#' scoresExample <- getScores(pathMEDExampleData, geneSets = "tmod", 
#'                              method = "GSVA")
#' @export
getScores <- function(inputData,
    geneSets,
    method = "GSVA",
    labels = NULL,
    cores = 1,
    use.assay = 1,
    ...) {
    if (is(inputData, "data.frame")) {
        inputData <- as.matrix(inputData)
    }

    if (is(inputData, "ExpressionSet")) {
        inputData <- Biobase::exprs(inputData)
    }

    if (is(inputData, "SummarizedExperiment")) {
        inputData <- as.matrix(SummarizedExperiment::assay(inputData, 
                                                            use.assay))
    }

    if (is.data.frame(geneSets)) {
        if (method %in% c("GSVA", "ssGSEA", "singscore", "Z-score", "Plage",
                            "M-Scores")) {
            stop("data frame based network is not permitted for method
                    selected")
        } else {
            if (all(c("source", "target", "weight") %in% colnames(geneSets))) {
                if (!"mor" %in% colnames(geneSets)) {
                    geneSets$mor <- rep(1, nrow(geneSets))
                }
            } else {
                stop(
                    "If a data frame based network is included in geneSets,
                    columns",
                    "source, target and weight are required"
                )
            }
        }
    }
    else if (is(geneSets, "GeneSetCollection")) {
        geneSets <- .gsc_to_list(geneSets)
    }
    else {
        if (!is.list(geneSets)) {
            data_env <- new.env(parent = emptyenv())
            data("genesetsData", envir = data_env, package = "pathMED")
            genesetsData <- data_env[["genesetsData"]]
            if (geneSets %in% names(genesetsData)) {
                geneSets <- genesetsData[[geneSets]]
            } else {
                stop("geneSets must be a list of genesets, a data frame based
                        network or a database name: ",
                    paste(names(genesetsData), collapse = ", ")
                )
            }
        }
    }

    params <- list(...)
    if (method %in% c("GSVA", "ssGSEA", "Z-score", "Plage")) {
        if (method == "GSVA") {
            params <- params[names(params) %in% c(
                "minSize", "maxSize",
                "kcdfNoneMinSampleSize",
                "tau",
                "maxDiff", "absRanking",
                "sparse", "checkNA", "use"
            )]
            paramMatrix <- do.call(GSVA::gsvaParam, c(
                list(inputData, geneSets,
                    kcdf = "Gaussian"
                ),
                params
            ))
        } else if (method == "ssGSEA") {
            params <- params[names(params) %in% c(
                "minSize", "maxSize", "alpha",
                "normalize", "checkNA",
                "use"
            )]
            paramMatrix <- do.call(GSVA::ssgseaParam, c(
                list(
                    inputData,
                    geneSets
                ),
                params
            ))
        } else if (method == "Z-score") {
            params <- params[names(params) %in% c("minSize", "maxSize")]
            paramMatrix <- do.call(GSVA::zscoreParam, c(
                list(
                    inputData,
                    geneSets
                ),
                params
            ))
        } else {
            params <- params[names(params) %in% c("minSize", "maxSize")]
            paramMatrix <- do.call(GSVA::plageParam, c(
                list(
                    inputData,
                    geneSets
                ),
                params
            ))
        }
        res <- GSVA::gsva(paramMatrix,
            BPPARAM = BiocParallel::SnowParam(workers = cores)
        )
        attr(res, "geneSets") <- NULL
    } else if (method == "singscore") {
        params <- params[names(params) %in% c(
            "subSamples", "centerScore",
            "dispersionFun", "knownDirection"
        )]
        rankMatrix <- singscore::rankGenes(inputData, tiesMethod = "average")
        listSing <- lapply(geneSets, function(x) {
            do.call(singscore::simpleScore, c(
                list(
                    rankData = rankMatrix,
                    upSet = x
                ),
                params
            ))
        })
        listScores <- lapply(listSing, function(x) x$TotalScore)

        if (is(listScores, "list")) {
            res <- do.call(rbind, listScores)
        } else {
            res <- t(listScores)
        }
        colnames(res) <- colnames(inputData)
    } else if (method %in% c(
        "AUCell", "MDT", "MLM", "ORA", "UDT", "ULM",
        "FGSEA", "norm_FGSEA", "WMEAN", "norm_WMEAN",
        "corr_WMEAN", "WSUM", "norm_WSUM", "corr_WSUM"
    )) {
        # Network necesary for these methods
        if (!is.data.frame(geneSets)) {
            net <- do.call("rbind", lapply(seq_len(length(geneSets)),
                                            function(x) {
                res <- data.frame(
                    "source" = rep(
                        names(geneSets)[[x]],
                        length(geneSets[[x]])
                    ),
                    "target" = as.character(geneSets[[x]]),
                    "weight" = rep(1, length(geneSets[[x]])),
                    "mor" = rep(1, length(geneSets[[x]]))
                )
                return(res)
            }))
        } else {
            net <- geneSets
        }

        if (method %in% c("AUCell", "MDT", "MLM", "ORA", "UDT", "ULM")) {
            if (method == "AUCell") {
                if (!requireNamespace("AUCell")) {
                    message("Install the AUCell package to run the AUCell
                            method")
                }
                params <- params[names(params) %in% c(
                    "aucMaxRank", "seed",
                    "minsize"
                )]
                scoreMatrix <- do.call(decoupleR::run_aucell, c(
                    list(
                        mat = inputData, network = net, .source = "source",
                        .target = "target"
                    ),
                    params
                ))
            } else if (method == "MDT") {
                params <- params[names(params) %in% c(
                    "center", "na.rm",
                    "trees", "min_n", "seed",
                    "minsize"
                )]
                scoreMatrix <- do.call(decoupleR::run_mdt, c(
                    list(
                        mat = inputData,
                        network = net,
                        nproc = cores,
                        .source = "source",
                        .target = "target",
                        .mor = "mor"
                    ),
                    params
                ))
            } else if (method == "MLM") {
                message("Filtering collinear pathways...")
                co.lin <- as.data.frame(decoupleR::check_corr(net,
                    .source = "source",
                    .target = "target",
                    .mor = "mor"
                ))
                net <- net[!net$source %in% as.character(co.lin[
                    co.lin$correlation > 0.7, "source"
                ]), ]

                params <- params[names(params) %in% c(
                    "center", "na.rm",
                    "minsize"
                )]
                scoreMatrix <- do.call(decoupleR::run_mlm, c(
                    list(
                        mat = inputData,
                        network = net,
                        .source = "source",
                        .target = "target",
                        .mor = "mor"
                    ),
                    params
                ))
            } else if (method == "ORA") {
                params <- params[names(params) %in% c(
                    "n_up", "n_bottom",
                    "n_background",
                    "with_ties",
                    "seed", "minsize"
                )]
                scoreMatrix <- do.call(decoupleR::run_ora, c(
                    list(
                        mat = inputData,
                        network = net,
                        .source = "source",
                        .target = "target"
                    ),
                    params
                ))
            } else if (method == "UDT") {
                params <- params[names(params) %in% c(
                    "center", "na.rm",
                    "min_n", "seed",
                    "minsize"
                )]
                scoreMatrix <- do.call(decoupleR::run_udt, c(
                    list(
                        mat = inputData,
                        network = net,
                        .source = "source",
                        .target = "target",
                        .mor = "mor"
                    ),
                    params
                ))
            } else if (method == "ULM") {
                params <- params[names(params) %in% c(
                    "center", "na.rm",
                    "minsize"
                )]
                scoreMatrix <- do.call(decoupleR::run_ulm, c(
                    list(
                        mat = inputData,
                        network = net,
                        .source = "source",
                        .target = "target",
                        .mor = "mor"
                    ),
                    params
                ))
            }

            scoreMatrix <- as.data.frame(scoreMatrix[, c(
                "source", "condition",
                "score"
            )])
        } else if (method %in% c("FGSEA", "norm_FGSEA")) {
            if (!requireNamespace("fgsea")) {
                message("Install the fgsea package to run the FGSEA or the 
                            norm_FGSEA method")
            }
            params <- params[names(params) %in% c("seed", "minsize", "times")]
            scoreMatrix <- do.call(decoupleR::run_fgsea, c(
                list(
                    mat = inputData,
                    network = net,
                    .source = "source",
                    .target = "target",
                    nproc = cores
                ),
                params
            ))

            if (method == "FGSEA") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "fgsea",
                        c("source", "condition", "score")
                    ]
                )
            } else {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "norm_fgsea",
                        c("source", "condition", "score")
                    ]
                )
            }
        } else if (method %in% c("WMEAN", "norm_WMEAN", "corr_WMEAN")) {
            params <- params[names(params) %in% c(
                "seed", "minsize", "sparse",
                "randomize_type", "times"
            )]
            scoreMatrix <- do.call(decoupleR::run_wmean, c(
                list(
                    mat = inputData,
                    network = net,
                    .source = "source",
                    .target = "target",
                    .mor = "mor"
                ),
                params
            ))

            if (method == "WMEAN") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "wmean",
                        c("source", "condition", "score")
                    ]
                )
            } else if (method == "norm_WMEAN") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "norm_wmean",
                        c("source", "condition", "score")
                    ]
                )
            } else {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "corr_wmean",
                        c("source", "condition", "score")
                    ]
                )
            }
        } else if (method %in% c("WSUM", "norm_WSUM", "corr_WSUM")) {
            params <- params[names(params) %in% c(
                "seed", "minsize", "sparse",
                "randomize_type", "times"
            )]
            scoreMatrix <- do.call(decoupleR::run_wsum, c(
                list(
                    mat = inputData,
                    network = net,
                    .source = "source",
                    .target = "target",
                    .mor = "mor"
                ),
                params
            ))

            if (method == "WSUM") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "wsum",
                        c("source", "condition", "score")
                    ]
                )
            } else if (method == "norm_WSUM") {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "norm_wsum",
                        c("source", "condition", "score")
                    ]
                )
            } else {
                scoreMatrix <- as.data.frame(
                    scoreMatrix[
                        scoreMatrix$statistic == "corr_wsum",
                        c("source", "condition", "score")
                    ]
                )
            }
        }

        scoreMatrix <- reshape2::dcast(scoreMatrix, source ~ condition,
            value.var = "score"
        )
        rownames(scoreMatrix) <- scoreMatrix$source
        scoreMatrix <- scoreMatrix[, -1]
        res <- scoreMatrix
    } else if (method == "M-Scores") {
        if (is.null(labels)) {
            stop("Labels parameter must be used for M-Scores method")
        }

        if (0 %in% labels) {
            HealthyData <- inputData[, labels == 0]
            PatientData <- inputData[, labels != 0]
            PatientData <- cbind(PatientData, HealthyData)
        } else if ("Healthy" %in% labels) {
            HealthyData <- inputData[, labels == "Healthy"]
            PatientData <- inputData[, labels != "Healthy"]
            PatientData <- cbind(PatientData, HealthyData)
        } else {
            stop("Reference samples in labels must be specified with 0
                or 'Healthy'")
        }

        message(
            "Healthy samples supplied. Calculating M-Scores using ",
            "healthy samples as reference for ", ncol(PatientData),
            " samples"
        )

        HealthyMean <- rowMeans(HealthyData, na.rm = TRUE)
        HealthySD <- apply(HealthyData, 1, function(x) {
            sd(x,
                na.rm = TRUE
            )
        })
        HealthyMeanSD <- cbind(HealthyMean, HealthySD)

        HealthyMeanSD <- HealthyMeanSD[HealthyMeanSD[, 2] != 0, ]

        res <- BiocParallel::bplapply(seq_len(ncol(PatientData)),
            function(column, geneSets,
    HealthyMeanSD) {
                Patient <- PatientData[, column,
                    drop = TRUE
                ]
                res.i <- lapply(geneSets,
                    .getMscorePath,
                    HealthyMeanSD =
                        HealthyMeanSD,
                    Patient = Patient
                )
                res.i <- as.data.frame(do.call(
                    "rbind", res.i
                ))
                return(res.i)
            },
            geneSets = geneSets,
            HealthyMeanSD = HealthyMeanSD,
            BPPARAM = BiocParallel::SnowParam(
                workers = cores, progressbar = TRUE
            )
        )

        res <- do.call("cbind", res)
        colnames(res) <- colnames(PatientData)
    }

    return(as.matrix(res))
}
