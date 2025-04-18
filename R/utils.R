## Function to retrieve Mscore for one patient and one path against healthy
## control distribution
.getMscorePath <- function(path,
    Patient,
    HealthyMeanSD) {
    genes.path <- Reduce(intersect, list(
        rownames(HealthyMeanSD),
        names(Patient),
        as.character(unlist(path))
    ))
    path.mscore <- 0 # default

    if (length(genes.path) > 2) {
        # at least 3 genes in each path
        tmpRef <- HealthyMeanSD[genes.path, ]
        tmpPat <- Patient[genes.path]

        # Zscore by path
        Zscore.genes <- (tmpPat - tmpRef[, 1]) / tmpRef[, 2]
        path.mscore <- mean(Zscore.genes, na.rm = TRUE)
    }
    names(path.mscore) <- names(path)
    return(path.mscore)
}



## Impute M-score for a patient from gene expression based on a gene expression
## reference of patients
# @ patient: vector of gene expression of a patient (names of vector must be
## gene symbols)
# @ Reference.normalized: gene expression matrix of patients normalized by
## z-score (by patient)
# @ Reference.mscore: matrix with M-scores of patients contained in
## Reference.normalized
# @ k: number of neighbours to impute mscore for patient
.getNearSample <- function(patient,
    Ref.norm,
    Ref.mscore,
    k = 5) {
    patient <- .normSamples(patient[intersect(
        names(patient),
        rownames(Ref.norm)
    )])

    distances <- apply(Ref.norm, 2, function(x) {
        stats::dist(rbind(patient, x))
    })
    names(distances) <- colnames(Ref.norm)

    distances <- distances[order(distances, decreasing = FALSE)]
    tmp.mscore <- Ref.mscore[, names(distances)[seq_len(k)]]
    tmp.mscore <- apply(tmp.mscore, 1, mean)

    return(list(
        "mscores" = tmp.mscore,
        "distance" = mean(distances[seq_len(k)], na.rm = TRUE)
    ))
}


## Normalize quantitative values (expression, M-scores) of a patient by Z-score
# @ x: Numeric vector of expression from a sample
.normSamples <- function(x) {
    x <- (x - (median(x, na.rm = TRUE))) / stats::sd(x, na.rm = TRUE)
    return(x)
}


## Function to build a normalized reference of patients with which to impute
## the m-scores for new samples
# @ expr.list: list with expression data from one or multiple studies (used as
## reference)
# @ geneset.list: output of DiseasePaths function (Pathways relevant for the
## disease)
# @ mscore.list: list with m-score matrices from the same studies that expr.list
.MReference <- function(expr.list,
                        geneset.list,
                        mscore.list) {
    ## Get common genes for all samples of the reference
    all.genes <-
        lapply(expr.list, function(i) {
            x <- rownames(i[[1]])
        })
    all.genes$geneset <- unique(as.character(unlist(geneset.list)))
    common.genes <- Reduce(intersect, all.genes)

    ## Reference gene-expression
    names(expr.list) <- NULL
    Reference <-
        lapply(expr.list, function(x) {
            return(x[[1]][common.genes, ])
        })
    Reference <- do.call(cbind, Reference)
    Reference.normalized <- apply(Reference, 2, .normSamples)

    ## Reference M-score
    Reference.mscore <- lapply(mscore.list, function(x) {
        return(x[names(geneset.list), ])
    })
    Reference.mscore <- do.call(cbind, Reference.mscore)

    return(
        list(
            Reference.mscore = Reference.mscore,
            Reference.normalized = Reference.normalized
        )
    )
}


## Function to remove messages and cat text from a function output
.removeOutText <- function(...) {
    tmpf <- tempfile()
    sink(tmpf)
    on.exit({
        sink()
        file.remove(tmpf)
    })
    out <- suppressMessages(eval(...))
    out
}


## Function to created class-balanced fold
.makeClassBalancedFolds <- function(y, kfold, repeats, varType, paired) {
    ## kfolds by samples (only one sample by patient) ··························
    if (is.null(paired)) {
        if (varType == "character") {
            y <- data.frame("value" = y, index = seq_len(length(y)))
            y <- split(y, y$value)

            splittedFolds <- lapply(y, function(x) {
                folds <- caret::createMultiFolds(x$index, k = kfold, repeats)
                res <- lapply(folds, function(i) {
                    x[i, ]$index
                })
            })
            listFolds <- lapply(
                seq_len(length(splittedFolds[[1]])),
                function(k) {
                    res <- c(unlist(lapply(
                        seq_len(length(splittedFolds)),
                        function(r) {
                            as.integer(unname(splittedFolds[[r]][k])[[1]])
                        }
                    )))
                }
            )
        } else {
            y <- data.frame("value" = y, index = seq_len(length(y)))
            folds <- caret::createMultiFolds(y$index, k = kfold, repeats)
            listFolds <- lapply(folds, function(i) {
                y[i, ]$index
            })
        }
        return(listFolds)

        ## kfolds by patients (several samples for each patient) ···············
    } else {
        if (varType == "character") {
            groups <- unique(y)
            ids <- lapply(groups, function(x) {
                return(unique(paired[y == x]))
            })
            names(ids) <- groups

            subsets <- lapply(seq_len(repeats), function(rp) {
                ## Random initiation
                tmp.ids <- lapply(ids, sample)
                tmp.ids <- do.call("rbind", lapply(tmp.ids, function(gr) {
                    return(data.frame(
                        "index" = gr,
                        "fold" = cut(seq_along(gr),
                            breaks = kfold,
                            labels = FALSE
                        )
                    ))
                }))
                rownames(tmp.ids) <- NULL

                ## Get folds of patients
                res <- list()
                for (f in unique(tmp.ids$fold)) {
                    filtered_indices <- tmp.ids$index[tmp.ids$fold != f]

                    sampleSel <- seq_len(length(paired))
                    res[[as.character(f)]] <- sampleSel[paired %in%
                        filtered_indices]
                }
                return(res)
            })
            subsets <- unlist(subsets, recursive = FALSE)
            names(subsets) <- NULL
            return(subsets)
        } else {
            ids <- unique(paired)

            subsets <- lapply(seq_len(repeats), function(rp) {
                gr <- ids[sample(seq_len(length(ids)), length(ids),
                    replace = FALSE
                )]
                tmp.ids <- data.frame(
                    "index" = gr,
                    "fold" = cut(seq_along(gr),
                        breaks = kfold,
                        labels = FALSE
                    )
                )

                res <- list()
                for (f in unique(tmp.ids$fold)) {
                    filtered_indices <- tmp.ids$index[tmp.ids$fold != f]

                    sampleSel <- seq_len(length(paired))
                    res[[as.character(f)]] <- sampleSel[paired %in%
                        filtered_indices]
                }
                return(res)
            })
            subsets <- unlist(subsets, recursive = FALSE)
            names(subsets) <- NULL

            return(subsets)
        }
    }
}


## Function to cluster pathways into co-expressed subpathways
.clusterPath <- function(data, path_name, minSplitSize, explainedVariance,
    maxSplits, cooccurrence = FALSE) {
    pca <- FactoMineR::PCA(t(data), graph = FALSE)
    pca_eig <- as.data.frame(pca$eig)
    pca_eig <-
        pca_eig[pca_eig$`cumulative percentage of variance` <
            explainedVariance, ]
    npcas <- nrow(pca_eig) + 1 ## Get K (npcas)
    if (all(is.na(pca_eig$`cumulative percentage of variance`))) {
        npcas <- 0
    }

    if (npcas > 1) {
        if (!is.null(maxSplits)) {
            if (npcas > maxSplits) {
                npcas <- maxSplits
            }
        }

        clust <- stats::kmeans(data, npcas)

        if (!all(table(clust$cluster) >= minSplitSize)) {
            ## Check if here are small clusters
            pca <- stats::prcomp(data, scale = FALSE, center = FALSE)

            ind <-
                factoextra::facto_summarize(pca,
                    element = "ind",
                    result = "coord",
                    axes = c(1, 2)
                )
            ind$cluster <- as.factor(clust$cluster)
            colnames(ind) <- c("name", "x", "y", "cood", "cluster")

            clusterPaths <- ind
            distance <- NULL

            # Distance between clusters (x,y)
            clusters <- as.numeric(unique(clusterPaths$cluster))
            clusters <- clusters[order(clusters, decreasing = FALSE)]
            for (cl in clusters) {
                distance <-
                    rbind(distance, c(
                        mean(clusterPaths$x[clusterPaths$cluster == cl],
                            na.rm = TRUE
                        ),
                        mean(clusterPaths$y[clusterPaths$cluster == cl],
                            na.rm = TRUE
                        )
                    ))
            }
            distance <-
                as.matrix(dist(
                    distance,
                    diag = NULL,
                    upper = TRUE,
                    method = "euclidean"
                ))
            diag(distance) <- NA
            rownames(distance) <- as.character(clusters)

            ## Joint small clusters
            for (iter in seq_len(minSplitSize)) {
                clusters <- table(clusterPaths$cluster)
                clusters <-
                    as.numeric(names(clusters[order(clusters,
                        decreasing = FALSE
                    )]))
                for (cl in clusters) {
                    if (as.numeric(table(clusterPaths$cluster == cl)["TRUE"]) <
                        minSplitSize) {
                        nearCl <- as.numeric(names(which.min(
                            distance[as.character(cl), ]
                        )))
                        clusterPaths[clusterPaths$cluster == cl, "cluster"] <-
                            nearCl
                        distance <-
                            distance[
                                !rownames(distance) %in% as.character(cl),
                                !colnames(distance) %in% as.character(cl)
                            ]
                    }
                }
                clusterPaths$cluster <-
                    factor(clusterPaths$cluster,
                        levels = unique(clusterPaths$cluster)
                    )
            }
        } else {
            clusterPaths <- data.frame(
                "name" = names(clust$cluster),
                "cluster" = as.numeric(clust$cluster)
            )
            rownames(clusterPaths) <- clusterPaths$name
            clusterPaths$cluster <-
                factor(clusterPaths$cluster,
                    levels = unique(clusterPaths$cluster)
                )
        }

        if (cooccurrence) {
            p.list <- data.frame("cluster" = clusterPaths$cluster)
            rownames(p.list) <- clusterPaths$name
            return(p.list)
        } else {
            paths <- split(clusterPaths, clusterPaths$cluster)
            p.list <- lapply(seq_len(length(paths)), function(pi) {
                genes <- as.character(rownames(paths[[pi]]))
            })
            names(p.list) <-
                paste0(path_name, ".split", seq(seq_len(length(p.list))))
            return(p.list)
        }
    } else {
        ## Only one K
        if (cooccurrence) {
            p.list <- data.frame("cluster" = rep(1, nrow(data)))
            rownames(p.list) <- rownames(data)
            return(p.list)
        } else {
            p.list <- list(rownames(data))
            names(p.list) <- path_name
            return(p.list)
        }
    }
}

# Transform GeneSetsCollection to list
.gsc_to_list <- function(gsc) {
    gene_sets <- lapply(gsc, GSEABase::geneIds)  # Extract gene IDs
    names(gene_sets) <- names(gsc)     # Assign names if available
    return(gene_sets)
}
