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
#' custom.tmod <- splitPathways(data=refData,
#' genesets="tmod",
#' customGeneset = NULL,
#' minPathSize = 10,
#' minSplitSize = 3,
#' maxSplits = 3,
#' explainedvariance = 50)
#' @export

splitPathways <- function (data, genesets, customGeneset = NULL, minPathSize = 10,
                           minSplitSize = 3, maxSplits = NULL, explainedvariance = 50)
{
    data.zscore <- lapply(data,function(x){
        Href <- data.frame("mean"=rowMeans(as.matrix(x$Healthy),na.rm = T),
                           "sd" = matrixStats::rowSds(as.matrix(x$Healthy),na.rm = T))

        x.zscore <- apply(as.matrix(x$Disease),2,function(pat){
            res <- (pat-Href$mean) / Href$sd
        })
        return(x.zscore)
    })
    ## Merge all samples
    genes <- Reduce(intersect,lapply(data.zscore,function(x){rownames(x)}))
    data.zscore <- lapply(data.zscore,function(x){x[genes,]})
    merged.datasets <- do.call("cbind",data.zscore)
    merged.datasets <- merged.datasets[,!colSums(is.na(merged.datasets)) > 0]

    ## Get pathway database
    if (genesets == "custom"){
        path.list <- customGeneset
    } else {
        path.list <- genesetsData[[genesets]]
    }

    new.path.list <- lapply(1:length(path.list),function(x){
        path_name <- names(path.list)[x]
        genes <- intersect(path.list[[x]], rownames(merged.datasets))
        if (length(genes) >= minPathSize) {
            tmp <- merged.datasets[genes,]
            pca <- FactoMineR::PCA(t(tmp),graph = F)

            pca_eig <- as.data.frame(pca$eig)
            pca_eig <- pca_eig[pca_eig$`cumulative percentage of variance` < explainedvariance,]
            npcas <- nrow(pca_eig) + 1

            if (length(npcas) != 0) {
                if (!is.null(maxSplits)) {
                    if (npcas > maxSplits) {
                        npcas = maxSplits
                    }
                }

                set.seed(123)
                x <- tmp
                clust <- stats::kmeans(x,npcas)
                data <- scale(x)
                cluster <- as.factor(clust$cluster)
                axes = c(1, 2)

                pca <- stats::prcomp(data, scale = FALSE, center = FALSE)
                ind <- factoextra::facto_summarize(pca, element = "ind", result = "coord",
                                                   axes = axes)
                ind$cluster <- cluster
                colnames(ind) <- c("name","x","y","cood","cluster")

                clusterPaths <- ind
                distance <- NULL
                clusters <- as.numeric(unique(clusterPaths$cluster))
                clusters <- clusters[order(clusters, decreasing = F)]
                for (d in 1:length(clusters)) {
                    distance <- rbind(distance, c(mean(clusterPaths$x[clusterPaths$cluster ==
                                                                          clusters[d]], na.rm = T), mean(clusterPaths$y[clusterPaths$cluster ==
                                                                                                                            clusters[d]], na.rm = T)))
                }
                distance <- as.matrix(dist(distance, diag = NULL,
                                           upper = T, method = "euclidean"))
                diag(distance) = NA
                rownames(distance) <- as.character(clusters)
                for (cl in 1:length(clusters)) {
                    if (as.numeric(table(clusterPaths$cluster ==
                                         clusters[cl])["TRUE"]) < minSplitSize) {
                        nearCl <- as.numeric(which.min(distance[as.character(clusters[cl]),
                        ]))
                        clusterPaths[clusterPaths$cluster == clusters[cl],
                                     "cluster"] <- nearCl
                        distance <- distance[!rownames(distance) %in%
                                                 as.character(clusters[cl]), ]
                    }
                }
                clusterPaths$cluster <- factor(clusterPaths$cluster)
                paths <- split(clusterPaths, clusterPaths$cluster)
                paths_list <- lapply(1:length(paths), function(i){
                    genes <- as.character(rownames(paths[[i]]))
                })
                names(paths_list) <- paste0(path_name,".split",seq(1:length(paths)))
                return(paths_list)
            } else{
                paths_list <- list(genes)
                names(paths_list) <- path_name
                return(paths_list)
            }
        } else{
            paths_list <- list(genes)
            names(paths_list) <- path_name
            return(paths_list)
        }
    })

    new_path_list <- list()
    for(x in 1:length(new.path.list)){
        if (length(new.path.list[[x]]) == 1){
            genes <- unlist(new.path.list[[x]])
            names(genes) <- NULL
            new_path_list[[names(new.path.list[[x]])]] <- genes
        } else{
            for (j in 1:length(new.path.list[[x]])){
                genes <- unlist(new.path.list[[x]][[j]])
                names(genes) <- NULL
                new_path_list[[names(new.path.list[[x]][j])]] <- genes
            }
        }
    }
    return(new_path_list)
}

