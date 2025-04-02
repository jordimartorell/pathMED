test_getScores <- function() {
    data(pathMEDExampleData)
    data(genesetsData)
    data(refData)
    reducedGeneSet <- genesetsData[["tmod"]][1:5]
    
    scoresExample <- getScores(refData[["dataset1"]], geneSets=reducedGeneSet,
                               labels= c(rep(1, 20), rep(0,10)),
                               method="M-Scores")
    checkEquals(round(scoresExample[1,1], 5), 6.02166)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="GSVA")
    checkEquals(round(scoresExample[1,1], 5), 0.60832)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="ssGSEA")
    checkEquals(round(scoresExample[1,1], 5), -0.28657)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="singscore")
    checkEquals(round(scoresExample[1,1], 5), -0.26854)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="Plage")
    checkEquals(round(scoresExample[1,1], 5), 0.16133)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="Z-score")
    checkEquals(round(scoresExample[1,1], 5), 1.54981)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="MDT")
    checkEquals(round(scoresExample[1,1], 5), 43.8289)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="MLM")
    checkEquals(round(scoresExample[1,1], 5), -2.18034)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="ORA")
    checkEquals(round(scoresExample[4,14], 5), 0.82403)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="UDT")
    checkEquals(round(scoresExample[1,1], 5), 0.77665)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="ULM")
    checkEquals(round(scoresExample[1,1], 5), -3.46543)
    scoresExample <- getScores(pathMEDExampleData[,1:5], geneSets=reducedGeneSet, 
                               method="FGSEA")
    checkEquals(round(scoresExample[1,1], 5), -0.36691)
    scoresExample <- getScores(pathMEDExampleData[,1:5], geneSets=reducedGeneSet, 
                               method="norm_FGSEA")
    checkEquals(round(scoresExample[1,1], 5), -1.44215)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="WMEAN")
    checkEquals(round(scoresExample[1,1], 5), 3.87817)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="norm_WMEAN")
    checkEquals(round(scoresExample[1,1], 5), -3.11012)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="corr_WMEAN")
    checkEquals(round(scoresExample[1,1], 5), 6.58889)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="WSUM")
    checkEquals(round(scoresExample[1,1], 5), 42.65983)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="norm_WSUM")
    checkEquals(round(scoresExample[1,1], 5), -3.11012)
    scoresExample <- getScores(pathMEDExampleData, geneSets=reducedGeneSet, 
                               method="corr_WSUM")
    checkEquals(round(scoresExample[1,1], 5), 72.47777)
    
}