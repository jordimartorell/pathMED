#' Establish reference data structure for input to the pathMED pipeline
#'
#' @param data A list of data frames or a single data frame with samples in
#'  columns and features in rows.
#' @param metadata A list of data frames or a single data frame with information
#'  for each sample. Samples in rows and variables in columns.
#' @param groupVar Character or list of characters indicating the column name of
#'  @metadata classifying the samples in control and cases. If several metadatas
#'   are provided a @groupVar can be specified for each metadata.
#' @param controlGroup Character or list of characters indicating which @groupVar
#'  level corresponds to the control group, usually healthy samples. All other
#'  samples will be considered as cases, usually disease samples. If
#'  several @groupVar are provided a @controlGroup can be specified
#'  for each @groupVar
#'
#' @return A list of reference data frames that serves as input
#'  for getMscoresRef and splitPathways functions.
#'
#' @author Iván Ellson, \email{ivan.ellson.l@@gmail.com }
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @seealso \code{\link{buildRefObject}}
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(reference_datasets)
#'
#' refData <- buildRefObject(
#'   data = list(dataset1, dataset2, dataset3, dataset4),
#'   metadata = list(metadata1, metadata2, metadata3, metadata4),
#'   groupVar = "group",
#'   controlGroup = "Healthy_sample")
#'
#' metadata <- rbind(metadata1, metadata2, metadata3, metadata4)
#' refData <- buildRefObject(
#'   data = list(dataset1, dataset2, dataset3, dataset4),
#'   metadata = metadata, # also works with a metadata for all datasets
#'   groupVar = "group",
#'   controlGroup = "Healthy_sample")
#'
#' @export
buildRefObject <- function(data, metadata, groupVar, controlGroup){
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
  if (!methods::is(controlGroup, "list")) {
    controlGroup <- list(controlGroup)
  }

  notLogDataset <- list()
  for (i in 1:length(data)) { # check log transformation
    qx <- as.numeric(quantile(data[[i]], c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
                              na.rm=TRUE))
    notLog <- (qx[5] > 100) ||(qx[6] - qx[1] > 50 &&
                                 qx[2] > 0) || (qx[2] > 0 &&
                                                  qx[2] < 1 &&
                                                  qx[4] > 1 && qx[4] < 2)
    notLogDataset <- append(notLogDataset, notLog)
    names(notLogDataset)[[i]] <-  paste0("dataset", i)
  }
  if (any(notLogDataset==TRUE)) {
    warning(paste0("The following expression datasets do not have a log2 distribution: ",
                   paste0(names(notLogDataset[notLogDataset==TRUE]),
                          collapse = ", "),
                   ". Please check that your data is normalized and transformed to log2."))
  }

  refData <- list()
  if (length(unique(c(length(data),length(metadata),length(groupVar),length(controlGroup))))==1) {
    for (i in 1:length(data)) {
      refData[[length(refData)+1]] <- list(
        Disease = as.data.frame(data[[i]][,intersect(colnames(data[[i]]), rownames(metadata[[i]])[!metadata[[i]][,groupVar[[i]]]==controlGroup[[i]]])]),
        Healthy = as.data.frame(data[[i]][,intersect(colnames(data[[i]]), rownames(metadata[[i]])[metadata[[i]][,groupVar[[i]]]==controlGroup[[i]]])]))
      names(refData)[i] <- paste0("dataset", i)
    }
  }
  if (length(metadata)==length(data) & length(metadata)>1 & length(groupVar)==length(controlGroup) & length(groupVar)==1) {
    for (i in 1:length(data)) {
      refData[[length(refData)+1]] <- list(
        Disease = as.data.frame(data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[i]])[!metadata[[i]][,groupVar[[1]]]==controlGroup[[1]]])]),
        Healthy = as.data.frame(data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[i]])[metadata[[i]][,groupVar[[1]]]==controlGroup[[1]]])]))
      names(refData)[i] <- paste0("dataset", i)
    }
  }
  if (length(metadata)!=length(data) & length(unique(c(length(metadata),length(groupVar),length(controlGroup),1)))==1) {
    for (i in 1:length(data)) {
      refData[[length(refData)+1]] <- list(
        Disease = as.data.frame(data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[1]])[!metadata[[1]][,groupVar[[1]]]==controlGroup[[1]]])]),
        Healthy = as.data.frame(data[[i]][intersect(colnames(data[[i]]),rownames(metadata[[1]])[metadata[[1]][,groupVar[[1]]]==controlGroup[[1]]])]))
      names(refData)[i] <- paste0("dataset", i)
    }
  }
  if (length(refData)==0) {
    stop("Input elements must be lists of datasets or individual datasets")
  }
  return(refData)
}
