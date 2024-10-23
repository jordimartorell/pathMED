test_diseasePaths <- function() {
    data(refData)
    refMscore <- getMscoresRef(data=refData,
                               genesets="tmod",
                               cores=1)
    checkEquals(refMscore$mscores[[1]][1,1], 6.021655503)
}
