#' Calculate M-scores for a dataset
#'
#' @param Patient Expression matrix of cases or numeric vector with one sample.
#' @param Healthy Expression matrix of healthy controls.
#' @param genesets Output from the diseasePaths function, or a list of pathways
#' @param nk Only for Healthy=NULL, Number of most similar samples to impute
#' M-scores.
#' @param maxDistance Maximum distance between patient's expression and k-samples
#' from patient-reference. Mscores for patients above this distance will not be imputed.
#' @param sqrtZmean Select if, when calculating the Mscores, divide by the 
#' square root of the number of genes when computing the mean of the Zscores. 
#' Boolean, FALSE by default.
#' @param cores Number of cores to be used.
#'
#' @return A list with the results of each of the analyzed regions. For each
#' region type, a data frame with the results and a list with the probes
#' associated to each region are generated. In addition, this list also contains
#' the input methData, pheno and platform objects
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{diseasePaths}}, \code{\link{getML}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(refData, exampleData)
#' exampleRefMScore <- getMscoresRef(data=refData, genesets="tmod")
#' relevantPaths <- diseasePaths(MRef=exampleRefMScore, min_datasets=3,
#' perc_samples=10)
#' MScoresExample <- getMscores(genesets = relevantPaths, Patient = exampleData,
#' Healthy = NULL, nk = 5)
#' @export
getMscores <- function(Patient,
                       Healthy=NULL,
                       genesets,
                       nk=5,
                       maxDistance= 30,
                       sqrtZmean=FALSE,
                       cores = 1){
  
  if (is.null(Healthy) & is.null(nk)) {
    stop("If Healthy=NULL, nk must be defined")
  }

  if(names(genesets)[2]=="reference"){
    path.list <- genesets[[1]]
    Reference <- genesets[[2]]
  }else{
    path.list <- genesets
    if(is.null(Healthy)){
      stop("No Healthy controls nor Reference included")
    }
  }
  
  
  if(is.vector(Patient)){ # Only one patient
    message("Calculating Mscores for one sample")
    
    if(is.null(Healthy)){
      message("No healthy samples supplied, calculating M-scores by ",
              "similarity with samples from reference")
      
      genes <- intersect(names(Patient),
                         rownames(Reference$Reference.normalized))
      Patient <- Patient[genes]
      
      res.l <- .getNearSample(patient=Patient,
                              Ref.norm=Reference$Reference.normalized[
                                genes,],
                              Ref.mscore=Reference$Reference.mscore,
                              k=nk)
      res <- as.data.frame(res.l$mscores)
      colnames(res)<-"Mscores"
      
      if(res.l$distance > maxDistance){
        message("Distance between patient's expression and k-samples from ",
                "Patient-reference are higher than maxDistance. Mscores ",
                "for this patients will not be imputed...")
        res<-NULL    
      }
      
    } else {
      message("Healthy samples supplied, calculating M-scores using ",
              "healthy samples as reference")
      H <- data.frame(apply(Healthy,1,function(x){mean(x,na.rm = T)}),
                      apply(Healthy,1,function(x){sd(x,na.rm = T)}))
      rownames(H)<-rownames(Healthy)
      H<-H[ifelse(H[,2]==0,F,T),]      
      res <- lapply(path.list, function(x) {
        .getMscorePath(x, Patient=Patient, Healthy=H, sqrtZmean=sqrtZmean)
      })
      res <- as.data.frame(do.call("rbind", res))
      colnames(res) <- "Mscores"
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
        names(x) <- genes
        res.l <- .getNearSample(patient=x,
                                Ref.norm=Reference$Reference.normalized[genes,],
                                Ref.mscore=Reference$Reference.mscore,
                                k=nk)
      })
      
      Mscores <- do.call("cbind",lapply(Mscores,function(m.x){
        if(m.x$distance <= maxDistance){
          return(m.x$mscores)    
        }
      }))
      
      if(!is.null(Mscores)){
        if(ncol(Patient) != ncol(Mscores)){
         message("Distance between expression of ", 
                 abs(ncol(Patient) - ncol(Mscores))," patients and k-samples ",
                 "from Patient-reference are higher than maxDistance. Mscores ",
                 "for these patients will not be imputed...")
        } 
      }else{
       message("Distance between expression of ",ncol(Patient)," patients and ",
       "k-samples from Patient-reference are higher than maxDistance. Mscores ",
       "for these patients will not be imputed...")
      }

      res <- Mscores
      
    } else {
      message("Healthy samples supplied. Calculating M-scores using ",
              "healthy samples as reference for ", ncol(Patient),
              " patients")
      H <- data.frame(apply(Healthy,1,function(x){mean(x,na.rm = T)}),
                      apply(Healthy,1,function(x){sd(x,na.rm = T)}))
      rownames(H)<-rownames(Healthy)
      H<-H[ifelse(H[,2]==0,F,T),]
      res <- BiocParallel::bplapply(seq_len(ncol(Patient)), function(column, geneNames,
                                                                     path.list, Healthy) {
        pat <- Patient[,column, drop=TRUE]
        names(pat) <- geneNames
        res.i <- lapply(path.list,
                        .getMscorePath,
                        Healthy=Healthy,
                        Patient=pat, 
                        sqrtZmean=sqrtZmean)
        res.i <- as.data.frame(do.call("rbind", res.i))
        return(res.i)
      },
      geneNames=rownames(Patient),
      path.list=path.list,
      Healthy=H,
      BPPARAM=BiocParallel::SnowParam(workers = cores, progressbar=TRUE))
      res <- do.call("cbind", res)
      colnames(res) <- colnames(Patient)
    }
  }
  return(res)
}

