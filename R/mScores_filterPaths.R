#' Find relevant pathways from the reference M-scores
#'
#' @param MRef output from the mScores_createReference function
#' @param min_datasets number of datasets that each pathway must meet the
#' perc_samples threshold
#' @param perc_samples minimun percentage of samples in a dataset in which a
#' pathway must be significant
#' @param Pcutoff P-value cutoff for significance
#' @param plotMetrics Plot number of significant pathways selected based on the different combination of perc_samples and min_datasets parameters
#'
#' @return A list with the selected pathways
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{mScores_createReference}}
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
#' exampleRefMScore <- mScores_createReference(data=refData, genesets="tmod")
#' relevantPaths <- mScores_filterPaths (MRef=exampleRefMScore, min_datasets=3,
#' perc_samples=10)
#' @export
mScores_filterPaths  <- function(MRef,
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
      plot2<-FALSE # For the moment, plot2 is not needed
      if (grepl("\\.split",names(MRef[["genesets"]])) & plot2) {  # any(endsWith(names(MRef[["genesets"]]), ".split1"))
        splitted_annotations.pre = names(MRef[["genesets"]])
        splitted_annotations.post = names(genesets)

        noSplitAnns.pre <- gsub(".split*.", "", splitted_annotations.pre)
        freqTable.pre <- as.data.frame(table(as.numeric(table(noSplitAnns.pre))))
        freqTable.pre$paths <- "Before selecion"
        freqTable.pre$Var1 <- paste0(as.character(freqTable.pre$Var1), " splits")
        freqTable.pre[nrow(freqTable.pre)+1,] <- c("All pathways (non-splitted)", sum(freqTable.pre$Freq), "Before selecion")
        freqTable.pre[nrow(freqTable.pre)+1,] <- c("All pathways (splitted)", length(splitted_annotations.pre), "Before selecion")
        freqTable.pre$Freq <- as.numeric(freqTable.pre$Freq)

        noSplitAnns.post  <- gsub(".split*.", "", splitted_annotations.post)
        freqTable.post  <- data.frame(Var1 = freqTable.pre$Var1, Freq = 0, paths = "After selection")
        freqTable.post$Freq[1:length(as.numeric(table(as.numeric(table(noSplitAnns.post)))))] <- as.numeric(table(as.numeric(table(noSplitAnns.post))))
        freqTable.post$Freq[nrow(freqTable.post)-1] <- sum(freqTable.post$Freq)
        freqTable.post$Freq[nrow(freqTable.post)] <- length(splitted_annotations.post)

        freqTable <- rbind(freqTable.pre, freqTable.post)
        colnames(freqTable) <- c("n.splits", "n.anns", "selection")
        freqTable$n.anns <- as.numeric(freqTable$n.anns)
        freqTable$n.splits[freqTable$n.splits=="1 splits"] <- "No splits"
        freqTable$selection <- factor(freqTable$selection, levels = c("Before selecion", "After selection"))
        freqTable$n.splits <- factor(freqTable$n.splits,
                                     levels = c("All pathways (splitted)", "All pathways (non-splitted)", "No splits",
                                                unique(freqTable$n.splits)[!unique(freqTable$n.splits) %in% c("All pathways (splitted)", "All pathways (non-splitted)", "No splits")]))

        P2 <- ggplot(freqTable, aes(x = n.splits, y = n.anns, fill = selection)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme_bw() +
          ylab("Number of pathways") +
          xlab("") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        Plots <- ggpubr::ggarrange(P1, NULL, P2, ncol = 1, align = "v", heights = c(3, 0.2, 2))
        plot(Plots)
      } else {plot(P1)}
    }
  }
  return(list(genesets=genesets, reference=reference))
}
