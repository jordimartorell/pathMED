#' Split pathways into coexpressed subpathways
#'
#' @param refData A refData object structure: a list of lists, each one with a
#' cases expression matrix and controls expression matrix
#' (named as Disease and Healthy). It can be constructed with the buildRefObject
#'  function. A list with one or more expression matrices, without controls, can
#'  also be used.
#' @param geneSets A named list with each
#' gene set or the name of one preloaded database (gobp, gomf, gocc,
#' kegg, reactome, tmod).
#' @param minSplitSize numeric, minimum number of genes in a subpathway.
#' Smaller splits will be merged with the closest coexpressed subpathway.
#' @param minPathSize numeric, minimum number of genes in a pathway to consider
#' splitting it.
#' @param maxSplits numeric, maximum number of subpathways for a pathway. If
#' NULL (default), there is not limit.
#' @param explainedVariance numeric, percentage of cumulative variance explained
#' within a pathway. This parameter is used to select the number of subdivisions
#' of a pathway that manage to explain at least the percentage of variance
#' defined by explainedVariance.
#' @param percSharedGenes numeric, minimum percentage of common genes across
#' datasets to merge them before clustering. If NULL or this percentage is not
#' reached, clustering is performed for each dataset independently and
#' consensus subpathways are obtained from co-occurrence across datasets.
#'
#' @return A list with the subpathways.
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{buildRefObject}}, \code{\link{mScores_createReference}},
#'  \code{\link{getScores}}
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
#' custom.tmod <- dissectDB(refData, geneSets="tmod")
#' @export

dissectDB <- function(refData,
                      geneSets,
                      minPathSize=10,
                      minSplitSize=3,
                      maxSplits=NULL,
                      explainedVariance=60,
                      percSharedGenes=90){

    ## 1. Get Z-scores by gene
    if (is(refData[[1]], "list")) {
        z.data <- lapply(refData, function(x) {
            Href <- data.frame("mean"=rowMeans(as.matrix(x$Healthy), na.rm=TRUE),
                               "sd"=matrixStats::rowSds(as.matrix(x$Healthy),
                                                        na.rm=TRUE))
            x.zscore <- as.data.frame(do.call("cbind",
                                              lapply(seq_len(ncol(x$Disease)),
                                                     function(pat) {
                                                         pat.i <- x$Disease[,pat]
                                                         names(pat.i) <- rownames(
                                                             x$Disease)
                                                         pat.i <- pat.i[
                                                             rownames(Href)]
                                                         return((pat.i - Href$mean)
                                                                / Href$sd)
                                                     })))
            colnames(x.zscore) <- colnames(x$Disease)
            x.zscore <- x.zscore[apply(x.zscore, 1, function(xi) {
                sum(is.na(xi))}) == 0,]
            return(x.zscore)
        })
    } else {
        z.data <- lapply(refData, function(x) {
            x.zscore <- t(scale(t(x)))
            x.zscore <- x.zscore[apply(x.zscore, 1, function(xi) {
                sum(is.na(xi))}) == 0,]
            return(x.zscore)
        })
    }


    ## 2. Getpathway database
    if (!is.list(geneSets)) {
        if(geneSets %in% names(genesetsData)){
            geneSets <- genesetsData[[geneSets]]
        }else{
            stop(paste("geneSets must be a list of genesets or a",
                       "database name: "),
                 paste(names(genesetsData), collapse=", "))
        }
    }

    ## Join all datasets if the specified percentage of genes are shared
    if(length(z.data) > 1 & !is.null(percSharedGenes)){
        exp.gr <- utils::combn(seq_len(length(z.data)), 2)
        sharedGenes <- unlist(lapply(seq_len(ncol(exp.gr)),
                                     function(it){
                                         x <- z.data[[exp.gr[1, it]]]
                                         y <- z.data[[exp.gr[2, it]]]
                                         return(min(c(length(intersect(
                                             rownames(x),
                                             rownames(y))
                                             ) / length(rownames(x)),
                                                      length(intersect(
                                                          rownames(x),
                                                          rownames(y))
                                                          ) / length(rownames(y)
                                                                     )))*100)
                                     }))

        if(all(sharedGenes >= percSharedGenes)){ ## Join all datasets
            genes.sd <- Reduce(intersect, lapply(z.data, function(x){
                rownames(x)}))
            genes.sd <- genes.sd[!is.na(genes.sd)]
            z.data <- lapply(z.data, function(x) {x[genes.sd,]})
            merged.datasets <- do.call("cbind", z.data)
            merged.datasets <- merged.datasets[,!colSums(is.na(merged.datasets))
                                               > 0]
            z.data <- list(merged.datasets)
            names(z.data)[1] <- "mergedData"
        }
    }

    ## 3. Dissect pathways
    message("This proccess can take time...")
    if(length(geneSets) > 1){
        pb=txtProgressBar(min=1, max=length(geneSets), initial=0, style=3)
    }
    new.geneSets <- list()
    for(p in seq_len(length(geneSets))) { ## Loop for each pathway
        path_name <- names(geneSets)[p]
        if(length(geneSets) > 1){setTxtProgressBar(pb,p)}
        if(length(z.data) == 1){ ## One dataset - kmeans clustering
            genes <- intersect(geneSets[[path_name]], rownames(z.data[[1]]))

            if(length(genes) > minPathSize){
                tmp <- z.data[[1]][genes,]
                p.list <- .clusterPath(data=tmp,
                                       path_name=path_name,
                                       minSplitSize=minSplitSize,
                                       explainedVariance=explainedVariance,
                                       maxSplits=maxSplits,
                                       cooccurrence=FALSE)
                new.geneSets <- c(new.geneSets,p.list)
            }
            else{
                ## Pathway with small size
                p.list <- list(genes)
                names(p.list) <- path_name
                #new.geneSets[count] <- p.list
                new.geneSets <- c(new.geneSets,p.list)
            }


        }
        else{ ## Multiple datasets - cooccurrence based kmeans clustering ······

            clusters.p <- lapply(seq_len(length(z.data)),
                                 function(d){ ## loop - datasets
                                     genes <- intersect(geneSets[[path_name]],
                                                        rownames(z.data[[d]]))

                                     if(length(genes) > minPathSize){
                                         tmp <- z.data[[d]][genes,]
                                         p.list <- .clusterPath(data=tmp,
                                                    path_name=path_name,
                                                    minSplitSize=minSplitSize,
                                                    explainedVariance=
                                                        explainedVariance,
                                                    maxSplits=maxSplits,
                                                    cooccurrence=TRUE)
                                     }
                                     else{
                                         ## Pathway with small size
                                         p.list <- data.frame("cluster"=
                                                                  rep(1, length(
                                                                      genes)))
                                         rownames(p.list) <- genes
                                         return(p.list)
                                     }
                                 })

            # Co-occurrence_matrix
            allgenes <- unique(unlist(lapply(clusters.p,
                                             function(cl){rownames(cl[1])})))

            if(length(allgenes) > minPathSize){

                cooccurrence_matrix <- matrix(0, nrow=length(allgenes),
                                              ncol=length(allgenes),
                                              dimnames=list(allgenes, allgenes))

                for(d in seq_along(clusters.p)){
                    cls <- unique(as.numeric(clusters.p[[1]][1][["cluster"]]))
                    for(cl in cls){
                        genes <- rownames(clusters.p[[d]])[
                            as.numeric(clusters.p[[d]]$cluster)==cl]
                        cooccurrence_matrix[genes,genes] <- cooccurrence_matrix[
                            genes,genes] + 1
                    }
                }

                p.list <- .clusterPath(data=cooccurrence_matrix,
                                       path_name=path_name,
                                       minSplitSize=minSplitSize,
                                       explainedVariance=explainedVariance,
                                       maxSplits=maxSplits,
                                       cooccurrence=FALSE)
                new.geneSets <- c(new.geneSets, p.list)
            }
            else{
                p.list <- list(allgenes)
                names(p.list) <- path_name
                new.geneSets <- c(new.geneSets, p.list)
            }
        }
    }
    if(length(geneSets) > 1) {
        close(pb)
    }

    return(new.geneSets)
}
