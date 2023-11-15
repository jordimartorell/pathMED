## Function to retrieve Mscore for one patient and one path against healthy
## control distribution
.getMscorePath <- function(path,
                           Patient,
                           Healthy){
    genes.path <- Reduce(intersect, list(rownames(Healthy),
                                         names(Patient),
                                         as.character(unlist(path))))
    path.mscore <- 0 # default

    if(length(genes.path) > 2){ # at least 3 genes in each path
        tmpRef <- Healthy[genes.path,]
        tmpPat <- Patient[genes.path]

        # Zscore by path
        Zscore.genes <- (tmpPat-tmpRef[,1])/tmpRef[,2]
        path.mscore <- mean(Zscore.genes, na.rm=TRUE)

    }
    names(path.mscore) <- names(path)
    return(path.mscore)
}



## Impute m-score for a patient from gene expression based on a gene expression
## reference of patients
#@ patient: vector of gene expression of a patient (names of vector must be gene
## symbols)
#@ Reference.normalized: gene expression matrix of patients normalized by
## z-score (by patient)
#@ Reference.mscore: matrix with mscores of patients contained in
## Reference.normalized
#@ k: number of neighbours to impute mscore for patient
.getNearSample <- function(patient,
                           Ref.norm,
                           Ref.mscore,
                           k=5){
    patient <- .normSamples(patient[intersect(names(patient),
                                              rownames(Ref.norm))])

    distances <- apply(Ref.norm, 2, function(x) {
        stats::dist(rbind(patient, x))
    })
    names(distances) <- colnames(Ref.norm)

    distances <- distances[order(distances,decreasing=FALSE)]
    tmp.mscore <- Ref.mscore[,names(distances)[seq_len(k)]]
    tmp.mscore <- apply(tmp.mscore,1,mean)

    return(list("mscores"=tmp.mscore,
                "distance"=mean(distances[seq_len(k)],na.rm=T)))
}


## Normalize quantitative values (expression, mscores) of a patient by z-score
#@ x: Numeric vector of expression from a sample
.normSamples <- function(x){
    x <- (x-(mean(x,na.rm=TRUE)))/stats::sd(x,na.rm=TRUE)
    return(x)
}


## Function to build a normalized reference of patients with which to impute
## the m-scores for new samples
#@ expr.list: list with expression data from one or multiple studies (used as
## reference)
#@ geneset.list: output of DiseasePaths function (Pathways relevant for the
## disease)
#@ mscore.list: list with m-score matrices from the same studies that expr.list
.MReference <- function(expr.list,
                        geneset.list,
                        mscore.list){

    ## Get common genes for all samples of the reference
    all.genes <- lapply(expr.list, function(i){x <- rownames(i[[1]])})
    all.genes$geneset<-unique(as.character(unlist(geneset.list)))
    common.genes <- Reduce(intersect, all.genes)

    ## Reference gene-expression
    names(expr.list) <- NULL
    Reference <- lapply(expr.list, function(x) {return(x[[1]][common.genes,])})
    Reference <- do.call(cbind, Reference)
    Reference.normalized<-apply(Reference, 2, .normSamples)

    ## Reference M-score
    Reference.mscore <- lapply(mscore.list, function(x) {
        return(x[names(geneset.list),])})
    Reference.mscore <- do.call(cbind, Reference.mscore)

    return(list(Reference.mscore=Reference.mscore,
                Reference.normalized=Reference.normalized))
}


## Function to remove messages and cat text from a function output
.removeOutText <- function(...){
    tmpf <- tempfile()
    sink(tmpf)
    on.exit({sink()
        file.remove(tmpf)})
    out <- suppressMessages(eval(...))
    out
}


## Function to created class-balanced fold
.makeClassBalancedFolds<-function(y,kfold,repeats,varType){
  
  if(varType=="character"){
    
    y<-data.frame("value"=y,index=1:length(y))
    y<-split(y, y$value)
    
    splittedFolds<-lapply(y,function(x){
      folds<-caret::createMultiFolds(x$index,k = kfold,repeats)
      res<-lapply(folds,function(i){
        x[i,]$index
      })
    })
    listFolds<-lapply(1:length(splittedFolds[[1]]),function(k){
      res<-c(unlist(lapply(1:length(splittedFolds),function(r){
        as.integer(unname(splittedFolds[[r]][k])[[1]])
      })))
    })
    
  }else{
    y<-data.frame("value"=y,index=1:length(y))
    folds<-caret::createMultiFolds(y$index,k = kfold,repeats)
    listFolds<-lapply(folds,function(i){
      y[i,]$index
    })
  }
  return(listFolds)
}


## Function to cluster pathways into co-expressed circuits

.clusterPath<-function(data,path_name,minSplitSize,explainedvariance,
                       maxSplits,cooccurrence =  FALSE){
  
  set.seed(1234)
  pca <- FactoMineR::PCA(t(data),graph = F)
  pca_eig <- as.data.frame(pca$eig)
  pca_eig <- pca_eig[pca_eig$`cumulative percentage of variance` < explainedvariance,]
  npcas <- nrow(pca_eig) + 1 ## Get K (npcas)
  
  if (npcas > 1) {
    if (!is.null(maxSplits)) {if (npcas > maxSplits) {npcas = maxSplits}}
    
    clust <- stats::kmeans(data,npcas)
    data.sc <- scale(data)
    pca <- stats::prcomp(data.sc, scale = FALSE, center = FALSE)
    ind <- factoextra::facto_summarize(pca, element = "ind", result = "coord",
                                       axes = c(1, 2))
    ind$cluster <- as.factor(clust$cluster)
    colnames(ind) <- c("name","x","y","cood","cluster")
    
    clusterPaths <- ind
    distance <- NULL
    
    # Distance between clusters (x,y)
    clusters <- as.numeric(unique(clusterPaths$cluster))
    clusters <- clusters[order(clusters, decreasing = F)]
    for (cl in clusters) {
      distance <- rbind(distance, c(mean(clusterPaths$x[clusterPaths$cluster == cl], na.rm = T),
                                    mean(clusterPaths$y[clusterPaths$cluster == cl], na.rm = T)))
    }
    distance <- as.matrix(dist(distance, diag = NULL,upper = T, method = "euclidean"))
    diag(distance) = NA
    rownames(distance) <- as.character(clusters)
    
    ## Joint small clusters
    for(iter in 1:minSplitSize){
      clusters <- table(clusterPaths$cluster)
      clusters <- as.numeric(names(clusters[order(clusters, decreasing = F)]))
      for (cl in clusters) {
        if (as.numeric(table(clusterPaths$cluster == cl)["TRUE"]) < minSplitSize) {
          nearCl <- as.numeric(names(which.min(distance[as.character(cl),])))
          clusterPaths[clusterPaths$cluster == cl,"cluster"] <- nearCl
          distance <- distance[!rownames(distance) %in% as.character(cl),
                               !colnames(distance) %in% as.character(cl)]
        }}
      clusterPaths$cluster <- factor(clusterPaths$cluster,levels=unique(clusterPaths$cluster))
    }
    
    if(cooccurrence){
      p.list<-data.frame("cluster"=clusterPaths$cluster)
      rownames(p.list)<-clusterPaths$name
      return(p.list)
      
    }else{
      paths <- split(clusterPaths, clusterPaths$cluster)
      p.list <- lapply(1:length(paths), function(pi){
        genes <- as.character(rownames(paths[[pi]]))})
      names(p.list) <- paste0(path_name,".split",seq(1:length(p.list)))
      return(p.list)
    }
    
  }else{ ## Only one K
    if(cooccurrence){
      p.list<-data.frame("cluster"=rep(1,nrow(data)))
      rownames(p.list)<-rownames(data)
      return(p.list)
      
    }else{
      p.list<-list(rownames(data))
      names(p.list)<-path_name
      return(p.list) 
    }
  }
}

