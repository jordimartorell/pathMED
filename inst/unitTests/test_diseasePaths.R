test_diseasePaths <- function() {
    data(refData)
    refMscore <- getMscoresRef(data=refData,
                               genesets="tmod",
                               cores=1)
    relevantPaths <- diseasePaths(MRef=refMscore, min_datasets=3,
                                  perc_samples=10)
    checkEquals(names(relevantPaths$genesets)[1], "LI.M4.2")
}
