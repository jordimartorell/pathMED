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
