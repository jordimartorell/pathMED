test_diseasePaths <- function() {
    data(exampleRefMScore)
    relevantPaths <- diseasePaths(MRef=exampleRefMScore, min_datasets=3,
                                  perc_samples=10)
    checkEquals(names(relevantPaths$genesets)[1], "LI.M4.2")
}
