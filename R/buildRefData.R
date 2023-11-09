#' Establish reference data structure for input to the pathMED pipeline
#' 
#' @param data A list of data frames or a single data frame with samples in 
#' columns and features in rows.
#' @param metadata A list of data frames or a single data frame with information
#'  for each sample. Samples in rows and variables in columns.
#' @param groupVar Character or list of characters indicating the column name of
#'  @metadata classifying the samples in healthy and disease. If several metadatas
#'   are provided a @groupVar can be specified for each metadata.
#' @param healthyGroup Character or list of characters indicating which @groupVar
#'  level corresponds to healthy samples. All other samples will be 
#'  considered as disease. If several @groupVar are provided a @healthyGroup
#'  can be specified for each @groupVar 
#' 
#' @return A list of reference Healthy and Disease data frames that serves as 
#' input for getMscoresRef and splitPathways functions.
#'
#' @author Daniel Toro-Dominguez, \email{daniel.toro@@genyo.es}
#' @author Jordi Martorell-Marugan, \email{jordi.martorell@@genyo.es}
#'
#' @seealso \code{\link{buildRefData}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#'
#' @export
buildRefData <- function(data, metadata, groupVar, healthyGroup){
  if (!methods::is(data, "list")) {
    message("Input data is not a list, treating it as a single dataset.")
    data <- list(data)
  }
  if (!methods::is(metadata, "list")) {
    metadata <- list(metadata)
  }
  if (!methods::is(groupVar, "list")) {
    groupVar <- list(groupVar)
  }
  if (!methods::is(healthyGroup, "list")) {
    healthyGroup <- list(healthyGroup)
  }
  
  refData <- list()
  if (length(unique(c(length(data),length(metadata),length(groupVar),length(healthyGroup))))==1) {
    for (i in 1:length(data)) {
      i = 1
      refData[[length(refData)+1]] <- list(
        Disease = data[[i]][,intersect(colnames(data[[i]]), rownames(metadata[[i]])[!metadata[[i]][,groupVar[[i]]]==healthyGroup[[i]]])],
        Healthy = data[[i]][,intersect(colnames(data[[i]]), rownames(metadata[[i]])[metadata[[i]][,groupVar[[i]]]==healthyGroup[[i]]])])
      names(refData)[i] <- paste0("dataset", i)
    }
  }
  if (length(metadata)==length(data) & length(metadata)>1 & length(groupVar)==length(healthyGroup) & length(groupVar)==1) {
    for (i in 1:length(data)) {
      refData[[length(refData)+1]] <- list(
        Disease = data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[i]])[!metadata[[i]][,groupVar[[1]]]==healthyGroup[[1]]])],
        Healthy = data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[i]])[metadata[[i]][,groupVar[[1]]]==healthyGroup[[1]]])])
      names(refData)[i] <- paste0("dataset", i)
    }
  }
  if (length(metadata)!=length(data) & length(unique(c(length(metadata),length(groupVar),length(healthyGroup),1)))==1) {
    for (i in 1:length(data)) {
      refData[[length(refData)+1]] <- list(
        Disease = data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[1]])[!metadata[[1]][,groupVar[[1]]]==healthyGroup[[1]]])],
        Healthy = data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[1]])[metadata[[1]][,groupVar[[1]]]==healthyGroup[[1]]])])
      names(refData)[i] <- paste0("dataset", i)
    }
  }
  if (length(refData)==0) {
    stop("Input elements must be lists of datasets or individual datasets")
  }
  return(refData)
}
