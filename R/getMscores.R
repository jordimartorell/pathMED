
## PATIENT Y HEALTHY AHORA SON PATHS, CAMBIA PARA QUE SEAN MATRICES

#' Calculate M-scores for a dataset
#'
#' Subtitle
#'
#' @param Patient Expression matrix of cases or numeric vector with one sample.
#' @param Healthy Expression matrix of healthy controls.
#' @param genesets Output from the diseasePaths function
#' @param nk Only for Healthy=NULL, Number of most similar samples to impute
#' M-scores.
#' @param cores Number of cores to be used.
#'
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated. In addition, this list also contains
#' the input methData, pheno and platform objects
#'
#' @author Daniel Toro-Domínguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
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
#' \dontrun{
#' DATA<-readRDS(file = paste0(getwd(),sep="","/data/datasets.rds"))
#' }
#' DATA.Mscore <- GetMscoresReferencet(genesets=genesets, data=DATA)
#' head(DATA.Mscore)
#' @export
GetMscores <- function(Patient,
                       Healthy=NULL,
                       genesets,
                       nk=5,
                       cores = 1){


    path.list <- genesets[[1]]
    Reference <- genesets[[2]]

    if(is.vector(Patient)){ # Only one patient
        message("Calculating Mscores for one sample")

        if(is.null(Healthy)){
            message("No healthy samples supplied, calculating M-scores by ",
                    "similarity with samples from reference")

            genes <- intersect(names(Patient),
                               rownames(Reference$Reference.normalized))
            Patient <- Patient[genes]

            res <- .getNearSample(patient=Patient,
                                Ref.norm=Reference$Reference.normalized[genes,],
                                Ref.mscore=Reference$Reference.mscore,
                                k=nk)
            res <- as.data.frame(res)
            colnames(res)<-"Mscores"

        } else {
            message("Healthy samples supplied, calculating M-scores using ",
                    "healthy samples as reference")
            res <- BiocParallel::bplapply(path.list, function(x) {
                .getMscorePath(x, Patient=Patient, Healthy=Healthy)
            }, BPPARAM=BiocParallel::MulticoreParam(workers = cores))
            res <- as.data.frame(do.call("rbind",res)); colnames(res)<-"Mscores"
        }

    } else { ## Several patients
        if (is.null(Healthy)) {
            message("No healthy samples supplied. Calculating M-scores by ",
                    "similarity with samples from reference for ",
                    ncol(Patient), " patients")

            genes <- intersect(rownames(Patient),
                             rownames(Reference$Reference.normalized))
            Patient <- Patient[genes,]

            Mscores <- pbapply::pbapply(Patient, 2, function(x) {
                .getNearSample(patient=x,
                               Ref.norm=Reference$Reference.normalized[genes,],
                               Ref.mscore=Reference$Reference.mscore,
                               k=nk)
            })

            res <- Mscores

        } else {
            message("Healthy samples supplied. Calculating M-scores using ",
                    "healthy samples as reference for ", ncol(Patient),
                    " patients")

            res <- pbapply::pbapply(Patient, 2, function(pat) {
                names(pat) <- rownames(Patient)
                res.i <- BiocParallel::bplapply(path.list,
                                                .getMscorePath,
                                                Healthy=Healthy,
                                                Patient=pat,
                                                BPPARAM=
                                                BiocParallel::MulticoreParam(
                                                    workers = cores))
                res.i <- as.data.frame(do.call("rbind", res.i))
                return(res.i)
            })
            res <- do.call("cbind", res)
            colnames(res) <- colnames(Patient)

            # res <- data.frame(matrix(nrow=length(path.list),
            #                          ncol=ncol(Patient)))
            # colnames(res) =  colnames(Patient)
            # rownames(res) <- names(path.list)
            #
            # res <- apply(Patient, 2, function(pat) {
            #     # utils::setTxtProgressBar(pb,i)
            #     names(pat) <- rownames(Patient)
            #     res.i <- BiocParallel::bplapply(.getMscorePath(path.list,
            #                                                  Healthy=Healthy,
            #                                                  Patient=pat),
            #                                   BPPARAM = MulticoreParam(workers = cores))
            #
            #     res.i <- as.data.frame(do.call("rbind", res.i))
            #     colnames(res.i) <- colnames(Patient)[i]
            #     return(res.i)
            # })


        }
    }
    return(res)
}
