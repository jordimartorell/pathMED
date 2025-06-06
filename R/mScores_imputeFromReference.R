#' Estimate M-scores for a dataset without healthy controls
#'
#' @param inputData Data matrix, data frame ExpressionSet or
#' SummarizedExperiment. Feature names must match the gene sets nomenclature. 
#' To use preloaded databases, they must be gene symbols.
#' @param geneSets A named list with each gene set,
#' or the name of one preloaded database (go_bp, go_cc, go_mf, kegg, reactome,
#' pharmgkb, lincs, ctd, disgenet, hpo, wikipathways, tmod)
#' or a GeneSetCollection.
#' @param externalReference External reference created with
#' the mScores_createReference function.
#' @param nk Number of
#' most similar samples from the external reference to impute M-scores.
#' @param distance.threshold Only samples that do not
#' surpass the mean Euclidean distance of distance.threshold (by default = 30)
#' with the external reference are imputed. If NULL,impute all samples.
#' @param cores Number of cores to be used.
#' @param use.assay If SummarizedExperiments are used, the number of the assay 
#' to extract the data.
#' 
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated. In addition, this list also contains
#' the input methData, pheno and platform objects
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{mScores_filterPaths}}, \code{\link{trainModel}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(refData, pathMEDExampleData)
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
#' refMScores <- mScores_createReference(refObject,
#'     geneSets = "tmod", cores = 1
#' )
#'
#' exampleMScores <- mScores_imputeFromReference(pathMEDExampleData,
#'     geneSets = "tmod",
#'     externalReference = refMScores,
#'     distance.threshold = 50
#' )
#' @export
mScores_imputeFromReference <- function(inputData,
    geneSets,
    externalReference,
    nk = 5,
    distance.threshold = 30,
    cores = 1,
    use.assay = 1) {
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

    if (is(geneSets, "GeneSetCollection")) {
        geneSets <- .gsc_to_list(geneSets)
    }
    else if (!is(geneSets, "list")) {
        data_env <- new.env(parent = emptyenv())
        data("genesetsData", envir = data_env, package = "pathMED")
        genesetsData <- data_env[["genesetsData"]]
        geneSets <- genesetsData[[geneSets]]
    }

    Patient <- inputData

    message(
        "Calculating M-Scores by ",
        "similarity with samples from external reference for ",
        ncol(inputData), " patients"
    )

    ## Get expression and mscores matrices from reference
    ref.expression <- lapply(externalReference$input, function(x) x$Disease)

    sharedGenes <- intersect(
        Reduce(intersect, lapply(ref.expression, rownames)),
        as.character(unlist(geneSets))
    )

    ref.expression <- as.data.frame(do.call(
        "cbind",
        lapply(ref.expression, function(m) {
            m[sharedGenes, , drop = FALSE]
        })
    ))
    ## normalization by patient to minimize impact of magnitude
    ## differences across sets
    ref.expression <- apply(ref.expression, 2, .normSamples)

    ref.mscores <- as.data.frame(do.call(
        "cbind",
        lapply(
            externalReference$mscores,
            function(m) {
                m[names(geneSets), , drop = FALSE]
            }
        )
    ))

    MscoresList <- pbapply::pbapply(Patient, 2, function(x) {
        names(x) <- rownames(Patient)
        res.l <- .getNearSample(
            patient = x,
            Ref.norm = ref.expression,
            Ref.mscore = ref.mscores,
            k = nk
        )
    })

    Mscores <- do.call("cbind", lapply(MscoresList, function(m.x) {
        return(m.x$mscores)
    }))
    Distances <- c(unlist(lapply(MscoresList, function(m.x) {
        return(m.x$distance)
    })))

    if (!is.null(distance.threshold)) {
        rem.paths <- unname(ifelse(Distances > distance.threshold, TRUE, FALSE))

        if (sum(rem.paths) > 0) {
            message(
                "Distance between expression of ",
                sum(rem.paths), "/", length(Distances), " patients ",
                "and k-samples from Patient-reference are higher ",
                "than ", distance.threshold, ". Mscores ",
                "for these patients will not be imputed..."
            )

            if (sum(rem.paths) == length(Distances)) {
                Mscores <- NULL
            } else {
                Mscores <- Mscores[, !rem.paths, drop = FALSE]
            }
        }
    }

    return(list(
        "Mscores" = Mscores,
        "Distances" = data.frame(
            "ID" = names(Distances),
            "Distance" = as.numeric(Distances)
        )
    ))
}
