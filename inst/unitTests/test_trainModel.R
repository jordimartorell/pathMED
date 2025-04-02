test_trainModel <- function() {
    data(pathMEDExampleData, pathMEDExampleMetadata)
    set.seed(123)
    scoresExample <- getScores(pathMEDExampleData, geneSets="tmod", 
                                method="Z-score")
    modelsList <- methodsML("svmLinear", outcomeClass="character")
    trainedModel <- trainModel(inputData=scoresExample,
                               metadata=pathMEDExampleMetadata,
                               var2predict="Response",
                               models=modelsList,
                               Koutter = 2,
                               Kinner = 2,
                               repeatsCV = 1)
    checkEquals(round(trainedModel$stats[1,1], 5), -0.04167)
}
