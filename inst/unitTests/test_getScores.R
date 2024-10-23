test_diseasePaths <- function() {
    data(refData, exampleData)
    refMscore <- createReference(datasetsList=refData,
                               geneSets="tmod",
                               cores=1)
    relevantPaths <- diseasePaths(MRef=refMscore, min_datasets=3,
                                  perc_samples=10)
    MScoresExample <- getMscores(geneSets = relevantPaths,
                                 inputData = exampleData,
                                 externalReference = refMscore, nk = 5)
    checkEquals(round(MScoresExample[1,1], 5), -0.02066)
}
