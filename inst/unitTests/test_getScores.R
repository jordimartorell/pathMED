test_getScores <- function() {
    data(exampleData)
    scoresExample <- getScores(exampleData, geneSets="tmod", method="GSVA")
    checkEquals(round(scoresExample[1,1], 5), 0.60832)
}
