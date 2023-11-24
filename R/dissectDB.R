#' Cluster pathways into pathway-components
#'
#' @param data A refData object structure: a list of lists, each one with a
#' cases expression matrix and controls expression matrix (named as Disease and Healthy)
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
#' defined by explainedvariance.
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
#' set.seed(1234)
#' custom.tmod <- dissectDB(data=refData,
#' genesets="tmod",
#' customGeneset = NULL,
#' minPathSize = 10,
#' minSplitSize = 3,
#' maxSplits = 3,
#' explainedvariance = 50)
#' @export

dissectDB<-function(data,genesets,customGeneset = NULL, minPathSize = 10,
                    minSplitSize = 3, maxSplits=NULL,explainedvariance = 60,
                    percSharedgenes = 90){
  
  ## 1. Get zscores by gene
  z.data<-lapply(data,function(x){
    Href <- data.frame("mean"=rowMeans(as.matrix(x$Healthy),na.rm = T),
                       "sd" = matrixStats::rowSds(as.matrix(x$Healthy),na.rm = T))
    x.zscore <- as.data.frame(do.call("cbind",lapply(1:ncol(x$Disease),function(pat){
      pat.i<-x$Disease[,pat]
      names(pat.i)<-rownames(x$Disease)
      pat.i<-pat.i[rownames(Href)]
      return((pat.i-Href$mean) / Href$sd)
    })))
    colnames(x.zscore)<-colnames(x$Disease)
    x.zscore<-x.zscore[apply(x.zscore,1,function(xi){sum(is.na(xi))})==0,]
    return(x.zscore)
  })
  
  ## 2. Getpathway database
  if(genesets == "custom"){
    path.list <- customGeneset
  }else{
    path.list <- genesetsData[[genesets]]
  }
  
  ## Joint all datasets if a certain percentage of genes are shared (i.e. 90%)
  if(length(z.data)>1 & !is.null(percSharedgenes)){
    exp.gr<-utils::combn(1:length(z.data),2)
    sharedGenes<-unlist(lapply(1:ncol(exp.gr),function(it){
      x<-z.data[[exp.gr[1,it]]]
      y<-z.data[[exp.gr[2,it]]]
      return(min(c(length(intersect(rownames(x),rownames(y)))/length(rownames(x)),
                   length(intersect(rownames(x),rownames(y)))/length(rownames(y))))*100)
    }))
    
    if(all(sharedGenes>=percSharedgenes)){ ## Joint all datasets
      genes.sd <- Reduce(intersect,lapply(z.data,function(x){rownames(x)}))
      genes.sd<-genes.sd[!is.na(genes.sd)]
      z.data <- lapply(z.data,function(x){x[genes.sd,]})
      merged.datasets <- do.call("cbind",z.data)
      merged.datasets <- merged.datasets[,!colSums(is.na(merged.datasets)) > 0]
      z.data<-list(merged.datasets)
      names(z.data)[1]<-"mergedData"
    }
  }
  
  ## 3. Disect pathways
  cat("This proccess can take time...\n")
  pb = txtProgressBar(min = 1, max = length(path.list), initial = 0)
  new.path.list<-list()
  for(p in 1:length(path.list)){ ## Loop for each pathway
    path_name <- names(path.list)[p]
    setTxtProgressBar(pb,p)
    if(length(z.data)==1){ ## One dataset - kmeans clustering ··················
      genes <- intersect(path.list[[path_name]], rownames(z.data[[1]]))
      
      if(length(genes) > minPathSize){
        tmp<-z.data[[1]][genes,]
        p.list<-.clusterPath(data=tmp,
                             path_name = path_name,
                             minSplitSize = minSplitSize,
                             explainedvariance = explainedvariance,
                             maxSplits = maxSplits,
                             cooccurrence = FALSE)
        new.path.list<-c(new.path.list,p.list)
      }else{
        ## Pathway with small size
        p.list <- list(genes)
        names(p.list)<-path_name
        #new.path.list[count]<-p.list
        new.path.list<-c(new.path.list,p.list)
      }
      
      
    }else{ ## Multiple datasets - ccooccurrence based kmeans clustering ········
      
      clusters.p<-lapply(1:length(z.data),function(d){ ## lapply loop - datasets
        genes <- intersect(path.list[[path_name]], rownames(z.data[[d]]))
        
        if(length(genes) > minPathSize){
          tmp<-z.data[[d]][genes,]
          p.list<-.clusterPath(data=tmp,
                               path_name = path_name,
                               minSplitSize = minSplitSize,
                               explainedvariance = explainedvariance,
                               maxSplits = maxSplits,
                               cooccurrence = TRUE)
        }else{
          ## Pathway with small size
          p.list<-data.frame("cluster"=rep(1,length(genes)))
          rownames( p.list)<-genes
          return(p.list)
        }
      }) ## lapply loop - datasets
      
      # Cooccurrence_matrix
      allgenes<-unique(unlist(lapply(clusters.p,function(cl){rownames(cl[1])})))
      
      if(length(allgenes)>minPathSize){
        
        cooccurrence_matrix <- matrix(0, nrow = length(allgenes), ncol = length(allgenes),
                                        dimnames = list(allgenes, allgenes))
          
          for(d in seq_along(clusters.p)){
            cls<-unique(as.numeric(clusters.p[[1]][1]$cluster))
            for(cl in cls){
              genes<-rownames(clusters.p[[d]])[as.numeric(clusters.p[[d]]$cluster)==cl]
              cooccurrence_matrix[genes,genes]<-cooccurrence_matrix[genes,genes]+1
            }
          }
        
        p.list<-.clusterPath(data=cooccurrence_matrix,
                             path_name = path_name,
                             minSplitSize = minSplitSize,
                             explainedvariance = explainedvariance,
                             maxSplits = maxSplits,
                             cooccurrence = FALSE)
        new.path.list<-c(new.path.list,p.list)
      }else{
        p.list <- list(allgenes)
        names(p.list)<-path_name
        #new.path.list[count]<-p.list
        new.path.list<-c(new.path.list,p.list)
      }
    }
  } ## Loop for each pathway
  
  return(new.path.list)
  
}
