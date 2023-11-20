#' Find relevant pathways from the reference M-scores
#'
#' @param MRef output from the getMscoresRef function
#' @param min_datasets number of datasets that each pathway must meet the
#' perc_samples threshold
#' @param perc_samples minimun percentage of samples in a dataset in which a
#' pathway must be significant
#' @param Pcutoff P-value cutoff for significance
#' @param plotMetrics Plot number of significant pathways selected based on the different combination of perc_samples and min_datasets parameters
#'
#' @return A list with the selected pathways
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{getMscoresRef}}
#'
#' @import ggplot2
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(refData)
#' exampleRefMScore <- getMscoresRef(data=refData, genesets="tmod")
#' relevantPaths <- diseasePaths(MRef=exampleRefMScore, min_datasets=3,
#' perc_samples=10)
#' @export
diseasePaths <- function(MRef,
                         min_datasets=round(length(MRef[[1]]) * 0.34),
                         perc_samples=10,
                         Pcutoff=0.05,
                         plotMetrics = TRUE){
  MScores <- MRef[[1]]
  genesets <- MRef[[2]]
  expr.list <- MRef[[3]]

  if (length(MRef[[1]])==1) {min_datasets <- 1}
  min_datasets <- min(min_datasets, length(MRef[[1]]))

  HighDys.perc <- lapply(MScores, function(dat) {
    apply(dat, 1, function(x) {
      values <- as.numeric(x[!is.na(x)])
      values <- (length(values[abs(values) >=
                                 abs(stats::qnorm(Pcutoff))])/
                   length(x))*100
      return(values)
    })
  })
  HighDys.perc <- do.call("cbind", HighDys.perc)

  selected.path <- apply(HighDys.perc, 1, function(x){
    nperc <- length(x[x > perc_samples])
    ntimes <- ifelse(nperc >= min_datasets, TRUE, FALSE)
    return(ntimes)
  })
  message("Selected paths: ", sum(selected.path))

  genesets <- genesets[selected.path]

  # Build the expression reference
  reference <- .MReference(expr.list=expr.list,
                           mscore.list=MScores,
                           geneset.list=genesets)

  ## plotMetrics
  if(plotMetrics){
   all<-expand.grid(1:length(MScores),1:100)
   colnames(all)<-c("min_datasets","perc_samples")

   nPaths<-NULL
   for(i in 1:nrow(all)){
     res <- apply(HighDys.perc, 1, function(x){
       nperc <- length(x[x > all$perc_samples[i]])
       ntimes <- ifelse(nperc >= all$min_datasets[i], TRUE, FALSE)
       return(ntimes)
     })
     res<- sum(res)
     nPaths<-c(nPaths,res)
   }
   all$selected_paths<-nPaths

   agg<-aggregate(selected_paths~perc_samples,all,sum)
   agg<-agg[agg$selected_paths>=5,] ## Minimal significant pathways across
   all<-all[all$perc_samples %in% unique(agg$perc_samples),]
   all$min_datasets<-paste("min_datasets",all$min_datasets)

    if(nrow(all)>0){
      P1<-ggplot(all,aes(x=perc_samples,y=selected_paths,group=min_datasets,color=min_datasets))+
       geom_line() + theme_bw()+
       geom_vline(xintercept = perc_samples, linetype = "dashed", color = "grey")+
       scale_x_continuous(breaks = round(seq(min(all$perc_samples), max(all$perc_samples), by = 5),1)) +
       scale_y_continuous(breaks = round(seq(min(all$selected_paths), max(all$selected_paths), by = as.integer(max(all$selected_paths)/10)),1))
      plot(P1)
    }
  }

  return(list(genesets=genesets, reference=reference))
}
