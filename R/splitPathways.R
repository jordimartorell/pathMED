#' Cluster pathways into pathway-components
#'
#' @param data list of lists, each one with a cases expression matrix and
#' controls expression matrix, always in this order.
#' @param genesets character, name of the preloaded database (gobp, gomf, gocc,
#' kegg, reactome, tmod) or 'custom' to provide the annotation.
#' @param customGeneset Only if genesets == 'custom'. A named list with each
#' gene set.
#' @param minSplitSize numeric, minimum number of genes in a pathway-split.
#' Smaller splits will join the closest split. 
#' @param minPathSize numeric, minimum number of genes in a pathway to consider 
#' split the pathway
#' @param maxSplits numeric, maximum number of subdivisions for a pathway. NULL, 
#' any restrictions
#' @param explainedvariance numeric, percentage of cumulative variance explained 
#' within a pathway. This parameter is used to select the number of subdivisions
#' of a pathway that manage to explain at least the percentage of variance 
#' defined by explainedvariance
#'
#' @return A list with all pathways and pathway-subdivisions and their involved
#' genes
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
#' data(refData)
#' custom.tmod<-splitPathway(data=refData,
#' genesets="tmod",
#' customGeneset = NULL,
#' minPathSize = 10,
#' minSplitSize = 3,
#' maxSplits = 3,
#' explainedvariance = 50)
#' @export

splitPathways<-function(data,
                      genesets,
                      customGeneset = NULL,
                      minPathSize = 10,
                      minSplitSize = 3,
                      maxSplits = NULL,
                      explainedvariance = 50){
  
  require("matrixStats")
  require("FactoMineR")
  require("factoextra")
  
  ## Expression to Z-score by gene
  data.zscore<-lapply(data,function(x){
    Href<-data.frame("mean"=rowMeans(as.matrix(x$Healthy),na.rm = T),
                     "sd" = matrixStats::rowSds(as.matrix(x$Healthy),na.rm = T))
    
    x.zscore <- apply(as.matrix(x$Disease),2,function(pat){
      res <- (pat-Href$mean) / Href$sd
    })
    return(x.zscore)
  })
  ## Merge all samples
  genes<-Reduce(intersect,lapply(data.zscore,function(x){rownames(x)}))
  data.zscore<-lapply(data.zscore,function(x){x[genes,]})
  merged.datasets<-do.call("cbind",data.zscore)
  merged.datasets<-merged.datasets[,!colSums(is.na(merged.datasets)) > 0]
  
  ## Get pathway database
  if (genesets == "custom"){
    path.list <- customGeneset
  } else {
    path.list <- genesetsData[[genesets]]
  }
  
  new.path.list<-list()
  npath<-1
  for(i in 1:length(path.list)){
    genes<-intersect(path.list[[i]],rownames(merged.datasets))
    
    if(length(genes) >= minPathSize){
      
      tmp<-merged.datasets[genes,]
      pca<-PCA(t(tmp),graph = F)
      
      ## Select number of components
      npcas <-nrow(pca$eig[ifelse(pca$eig[,3] < explainedvariance,T,F),])+1

      if(length(npcas)!=0){
        
        if(!is.null(maxSplits)){
          if(npcas>maxSplits){npcas=maxSplits}
        }
        
        clusterPaths<-eclust(tmp, "kmeans", k=npcas)$clust_plot$data
        invisible(dev.off())
        
        ## Join small clusters to nearest cluster
        distance<-NULL
        clusters<-as.numeric(unique(clusterPaths$cluster))
        clusters<-clusters[order(clusters,decreasing = F)]
        for(d in 1:length(clusters)){
          distance<-rbind(distance, c(mean(clusterPaths$x[clusterPaths$cluster==clusters[d]],na.rm = T),
                             mean(clusterPaths$y[clusterPaths$cluster==clusters[d]],na.rm = T)))
          
        }
        distance<-as.matrix(dist(distance,diag = NULL,upper = T,method = "euclidean"))
        diag(distance)=NA
        rownames(distance)<-as.character(clusters)
        
        for(cl in 1:length(clusters)){
          if(as.numeric(table(clusterPaths$cluster==clusters[cl])["TRUE"])<minSplitSize){
            nearCl<-as.numeric(which.min(distance[as.character(clusters[cl]),]))
            clusterPaths[clusterPaths$cluster==clusters[cl],"cluster"]<-nearCl
            distance<-distance[!rownames(distance) %in% as.character(clusters[cl]),]
          }
        }
        clusterPaths$cluster<-factor(clusterPaths$cluster)
        
        paths<-split(clusterPaths,clusterPaths$cluster)
        for(P in 1:length(paths)){
          new.path.list[[npath]]<-as.character(rownames(paths[[P]]))
          names(new.path.list)[npath]<-paste0(names(path.list)[i],sep=".split.",names(paths[P]))
          npath<-npath+1
        }
        
      }else{
        new.path.list[[npath]]<-path.list[[i]]
        names(new.path.list)[npath]<-names(path.list)[i]
        npath<-npath+1
      }
      
    }else{
      new.path.list[[npath]]<-path.list[[i]]
      names(new.path.list)[npath]<-names(path.list)[i]
      npath<-npath+1
    }
  }
  
  return(new.path.list)
}

