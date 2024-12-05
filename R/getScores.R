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

    if (!is.list(geneSets)) {
        if(geneSets %in% names(genesetsData)){
            geneSets<-genesetsData[[geneSets]]
        }else{
            stop(paste("geneSets must be a list of genesets or a database name:",
               names(genesetsData)))
        }
    } 

    params<-list(...)
    if (method %in% c("GSVA", "ssGSEA", "Z-score", "Plage")) {
        if (method == "GSVA") {
            params<-params[names(params) %in% c("minSize","maxSize","kcdfNoneMinSampleSize",
                                     "tau","maxDiff","absRanking","sparse","checkNA","use")]
            paramMatrix <-  do.call(GSVA::gsvaParam, c(list(inputData, geneSets, kcdf = "Gaussian"), params))
        }
        else if (method == "ssGSEA") {
            params<-params[names(params) %in% c("minSize", "maxSize", "alpha", "normalize",
                                                "checkNA","use")]
            paramMatrix <-  do.call(GSVA::ssgseaParam, c(list(inputData, geneSets), params))
        }
        else if (method == "Z-score") {
            params<-params[names(params) %in% c("minSize", "maxSize")]
            paramMatrix <-  do.call(GSVA::zscoreParam, c(list(inputData, geneSets), params))
        }
        else {
            params<-params[names(params) %in% c("minSize", "maxSize")]
            paramMatrix <-  do.call(GSVA::plageParam, c(list(inputData, geneSets), params))
        }
        res <- GSVA::gsva(paramMatrix,
                          BPPARAM=BiocParallel::SnowParam(workers=cores))
    }

    else if (method == "singscore") {
        params<-params[names(params) %in% c("subSamples","centerScore","dispersionFun","knownDirection")]
        rankMatrix <- singscore::rankGenes(inputData, tiesMethod = "average")
        listSing <- lapply(geneSets, function(x) {
            do.call(singscore::simpleScore, c(list(rankData = rankMatrix, upSet = x), params))
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
        # Filter pathways by collinearity
        #' @param collinearity_threshold must be a number between 0 and 1
        if("collinearity_threshold" %in% names(params)){
            cat("Filtering collinear pathways...")
            co.lin<-as.data.frame(decoupleR::check_corr(net,.source = "source",.target = "target",.mor = "mor",))
            net<-net[!net$source %in% as.character(co.lin[co.lin$correlation > 
                                                            as.numeric(params[["collinearity_threshold"]]),"source"]),]
        }    

        if (method %in% c("AUCell", "MDT", "MLM", "ORA", "UDT", "ULM")) {

            if (method == "AUCell") {
                params<-params[names(params) %in% c("aucMaxRank","seed","minsize")]
                scoreMatrix <-do.call(decoupleR::run_aucell, c(list(mat=inputData, network=net,
                                       .source="source",.target="target"),params))
            }

            else if (method == "MDT") {
                params<-params[names(params) %in% c("center","na.rm","trees","min_n","seed","minsize")]
                scoreMatrix <-do.call(decoupleR::run_mdt, c(list(mat=inputData, network=net,
                                       nproc=cores,.source="source",.target="target",.mor="mor"),params))
            }

            else if (method == "MLM") {
                params<-params[names(params) %in% c("center","na.rm","minsize")]
                scoreMatrix <-do.call(decoupleR::run_mlm, c(list(mat=inputData, network=net,
                                       .source="source",.target="target",.mor="mor"),params))
            }

            else if (method == "ORA") {
                params<-params[names(params) %in% c("n_up","n_bottom","n_background","with_ties","seed",
                                                    "minsize")]
                scoreMatrix <-do.call(decoupleR::run_ora, c(list(mat=inputData, network=net,
                                       .source="source",.target="target"),params))
            }

            else if (method == "UDT") {
                params<-params[names(params) %in% c("center","na.rm","min_n","seed","minsize")]
                scoreMatrix <-do.call(decoupleR::run_udt, c(list(mat=inputData, network=net,
                                       .source="source",.target="target",.mor="mor"),params))
            }

            else if (method == "ULM") {
                params<-params[names(params) %in% c("center","na.rm","minsize")]
                scoreMatrix <-do.call(decoupleR::run_ulm, c(list(mat=inputData, network=net,
                                      .source="source",.target="target",.mor="mor"),params))
            }

            scoreMatrix <- as.data.frame(scoreMatrix[,c("source","condition","score")])
          
        }

        else if (method %in% c("FGSEA", "norm_FGSEA")) {
            params<-params[names(params) %in% c("seed","minsize")]
            scoreMatrix <-do.call(decoupleR::run_fgsea, c(list(mat=inputData, network=net,
                                  .source="source",.target="target", nproc=cores, times = 1,),params))

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
            params<-params[names(params) %in% c("seed","minsize","sparse","randomize_type")]
            scoreMatrix <-do.call(decoupleR::run_wmean, c(list(mat=inputData, network=net,
                                  .source="source",.target="target",.mor="mor", nproc=cores, times = 1,),params))

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
            params<-params[names(params) %in% c("seed","minsize","sparse","randomize_type")]
            scoreMatrix <-do.call(decoupleR::run_wsum, c(list(mat=inputData, network=net,
                                  .source="source",.target="target",.mor="mor", times = 1,),params))

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

