test_diseasePaths <- function() {
    data(pathMED)
    relevantPaths <- diseasePaths(MRef=refMscore, min_datasets=3,
                                  perc_samples=10)
    checkEquals(names(relevantPaths$genesets)[1], "LI.M4.2")
}
