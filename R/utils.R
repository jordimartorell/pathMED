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
        Zscore.genes <- abs((tmpPat-tmpRef[,1])/tmpRef[,2])
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

    return(tmp.mscore)
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

# Prepare parameters for tuning the different methods
.methodsML <- function(algorithms,
                       add=NULL,
                       outcomeClass,
                       training){
    availableMethods <- list(numeric=c('glm', 'lm', 'xgbTree', 'rf', 'knn',
                                       'nnet', 'svmRadial', 'svmLinear', 'lars',
                                       'rpart','gamboost', 'brnn', 'enet'),
                             character=c('glm', 'lda', 'xgbTree', 'rf', 'knn',
                                         'nnet', 'svmLinear','svmRadial', 'nb',
                                         'ada', 'gamboost'))
    methodList <- list(
        lm=caretModelSpec(method='lm'),
        glm=caretModelSpec(method='glm'),
        lda=caretModelSpec(method='lda'),
        xgbTree=caretModelSpec(method='xgbTree',
                               tuneGrid=expand.grid(.max_depth = c(2, 3, 4, 5,
                                                                   6, 8, 10),
                                                    .nrounds = 50,
                                                    .eta = c(0.01, 0.05, 0.1),
                                                    .gamma = c(0,1),
                                                    .colsample_bytree=c(0.1,
                                                                        0.4),
                                                    .min_child_weight=c(1, 10),
                                                    .subsample=c(0.5, 1))),
        rf=caretModelSpec(method='rf',
                          tuneGrid=data.frame(.mtry=seq(2,
                                                        round(sqrt(ncol(
                                                            training)-1)) ,1))),
        knn=caretModelSpec(method='knn', tuneGrid=data.frame(.k=seq(2, 50, 2))),
        svmLinear=caretModelSpec(method='svmLinear', tuneLength = 15),
        svmRadial=caretModelSpec(method='svmRadial', tuneLength = 15),
        nnet=caretModelSpec(method='nnet',
                            tuneGrid=expand.grid(.size=seq(from=1, to=4, by=1),
                                                 .decay=seq(from=0.1, to=0.5,
                                                            by=0.1))),
        nb=caretModelSpec(method='nb',
                          tuneGrid=expand.grid(.fL=c(0, 0.5, 1.0),
                                               .usekernel=TRUE,
                                               .adjust=c(0.5, 1.0))),
        lars=caretModelSpec(method='lars',
                            tuneGrid=expand.grid(.fraction=seq(.01, .99,
                                                               length=40))),
        rpart=caretModelSpec(method='rpart',
                             tuneGrid=expand.grid(.cp=seq(0, .02, .0001))),
        ada=caretModelSpec(method='ada',
                           tuneGrid=expand.grid(.maxdepth=25,
                                                .nu=2, .iter=100)),
        gamboost=caretModelSpec(method='gamboost',
                                tuneGrid=expand.grid(.mstop=c(50, 100, 150, 200,
                                                              250, 300),
                                                     .prune=c('yes', 'no'))),
        brnn=caretModelSpec(method='brnn', tuneGrid=expand.grid(
            .neurons=c(2:16))),
        enet=caretModelSpec(method='enet',
                            tuneGrid=expand.grid(.lambda=c(0, 0.005, 0.009,
                                                           0.05, 0.1),
                                                 .fractions=c(0.05, 0.09,0.5)))
    )

    if('all' %in% algorithms){
        algorithms <- c('glm', 'lm', 'lda', 'xgbTree', 'rf', 'knn', 'svmLinear',
                        'svmRadial', 'nnet', 'nb', 'lars', 'rpart', 'ada',
                        'gamboost', 'brnn', 'enet')
    }
    algorithms <- algorithms[algorithms %in% unlist(
        availableMethods[outcomeClass])]
    methodList <- methodList[algorithms]
    methodList <- c(methodList, add)

    return(methodList)
}

## Function to remove co-linear features
#@ data: matrix with features in columns and samples in rows. Column of group
## is required (outcome to predict)
#@ thresh: SU threshold
.fast.cor.FS<-function(data,
                       thresh=0.0025){
    stopifnot('group' %in% colnames(data))

    y <- as.factor(data$group)
    x <- subset(data, select = -c(group))
    dis <- FCBF::discretize_exprs(t(x))

    fcbf.res <- FCBF::fcbf(dis, y, verbose=TRUE, thresh)
    xx <- x[,fcbf.res$index]
    xx <- as.data.frame(cbind(group=y, xx))

    return(xx)
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
