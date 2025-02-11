test_trainModel <- function() {
    data(exampleData, exampleMetadata)
    set.seed(123)
    scoresExample <- getScores(exampleData, geneSets="tmod", method="GSVA")
    modelsList <- methodsML("svmLinear", outcomeClass="character")
    trainedModel <- trainModel(inputData=scoresExample,
                               metadata=exampleMetadata,
                               var2predict="Response",
                               models=modelsList,
                               Koutter = 2,
                               Kinner = 2,
                               repeatsCV = 1)
    checkEquals(round(trainedModel$stats[1,1], 5), 0.15119)
}
