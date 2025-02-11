#' Estimate M-scores for a dataset without healthy controls
#'
#' @param inputData Data matrix or data frame.
#' @param geneSets A named list with each
#' gene set or the name of one preloaded database (gobp, gomf, gocc,
#' kegg, reactome, tmod).
#' @param externalReference External reference created with
#' the mScores_createReference function.
#' @param nk Number of
#' most similar samples from the external reference to impute M-scores.
#' @param imputeAll By default, only samples that do not
#' surpass the distance of 30 with the external reference are imputed. If TRUE,
#' impute all samples.
#' @param cores Number of cores to be used.
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
#' \donttest{
#' data(reference_datasets, exampleData)
#'
#' refData <- buildRefObject(
#'     data = list(dataset1, dataset2, dataset3, dataset4),
#'     metadata = list(metadata1, metadata2, metadata3, metadata4),
#'     groupVar = "group",
#'     controlGroup = "Healthy_sample"
#' )
#'
#' refMScores <- mScores_createReference(refData,
#'     geneSets = "tmod", cores = 1
#' )
#'
#' exampleMScores <- mScores_imputeFromReference(exampleData,
#'     geneSets = "tmod",
#'     externalReference = refData
#' )
#' }
#' @export
mScores_imputeFromReference <- function(inputData,
    geneSets,
    externalReference,
    nk = 25,
    imputeAll = FALSE,
    cores = 1) {
    if (is(inputData, "data.frame")) {
        inputData <- as.matrix(inputData)
    }

    if (!is(geneSets, "list")) {
        geneSets <- genesetsData[[geneSets]]
    }


    Patient <- inputData

    message(paste(
        "Calculating M-Scores by",
        "similarity with samples from external reference for",
        ncol(inputData), "patients"
    ))

    genes <- intersect(
        rownames(geneSets),
        rownames(externalReference$Reference.normalized)
    )
    geneSets <- geneSets[genes, ]

    Mscores <- pbapply::pbapply(geneSets, 2, function(x) {
        names(x) <- genes
        res.l <- .getNearSample(
            patient = x,
            Ref.norm = externalReference$
                Reference.normalized[genes, ],
            Ref.mscore = externalReference$
                Reference.mscore,
            k = nk
        )
    })

    Mscores <- do.call("cbind", lapply(Mscores, function(m.x) {
        if (imputeAll == TRUE) {
            return(m.x$mscores)
        } else if (m.x$distance <= 30) {
            return(m.x$mscores)
        }
    }))

    if (!is.null(Mscores)) {
        if (ncol(PatientData) != ncol(Mscores)) {
            message(paste(
                "Distance between expression of",
                abs(ncol(PatientData) - ncol(Mscores)), "patients",
                "and k-samples from Patient-reference are higher",
                "than 30. Mscores",
                "for these patients will not be imputed..."
            ))
        }
    } else {
        message(paste(
            "Distance between expression of", ncol(PatientData),
            "patients and",
            "k-samples from Patient-reference are higher",
            "than 30. Mscores",
            "for these patients will not be imputed..."
        ))
    }

    return(Mscores)
}
