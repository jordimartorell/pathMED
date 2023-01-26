test_diseasePaths <- function() {
    data(refData, exampleData)
    refMscore <- getMscoresRef(data=refData,
                               genesets="tmod",
                               cores=1)
    relevantPaths <- diseasePaths(MRef=refMscore, min_datasets=3,
                                  perc_samples=10)
    MScoresExample <- getMscores(genesets = relevantPaths,
                                 Patient = exampleData,
                                 Healthy = NULL, nk = 5)
    checkEquals(round(MScoresExample[1,1], 5), -0.02066)
}
