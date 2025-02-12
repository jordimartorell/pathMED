test_predictExternal <- function() {
    data(reference_datasets)
    commonGenes <- intersect(rownames(dataset1), rownames(dataset2))
    dataset1 <- dataset1[commonGenes, ]
    dataset2 <- dataset2[commonGenes, ]
    
    scoresExample <- getScores(dataset1, geneSets="tmod", method="Z-score")
   
    set.seed(123)
    trainedModel <- trainModel(inputData=scoresExample,
                                metadata=metadata1,
                                var2predict="group",
                                models=methodsML("svmLinear",
                                                 outcomeClass="character"),
                                Koutter=2,
                                Kinner=2,
                                repeatsCV=1)
                                
    externalScores <- getScores(dataset2, geneSets="tmod", method="Z-score")
    realValues <- metadata2$group
    names(realValues) <- rownames(metadata2)
    predictions <- predictExternal(externalScores, trainedModel, 
                                   realValues=realValues)
    
    checkEquals(round(predictions$stats[1,2], 5), 0.96667)
}
