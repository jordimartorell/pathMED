test_predictExternal <- function() {
    data(refData)
    commonGenes <- intersect(rownames(refData$dataset1), rownames(refData$dataset2))
    dataset1 <- refData$dataset1[commonGenes, ]
    dataset2 <- refData$dataset2[commonGenes, ]
    
    scoresExample <- getScores(dataset1, geneSets="tmod", method="Z-score")
   
    set.seed(123)
    trainedModel <- trainModel(inputData=scoresExample,
                                metadata=refData$metadata1,
                                var2predict="group",
                                models=methodsML("svmLinear",
                                                 outcomeClass="character"),
                                Koutter=2,
                                Kinner=2,
                                repeatsCV=1)
                                
    externalScores <- getScores(dataset2, geneSets="tmod", method="Z-score")
    realValues <- refData$metadata2$group
    names(realValues) <- rownames(refData$metadata2)
    predictions <- predictExternal(externalScores, trainedModel, 
                                   realValues=realValues)
    
    checkEquals(round(predictions$stats[1,2], 5), 0.96667)
}
