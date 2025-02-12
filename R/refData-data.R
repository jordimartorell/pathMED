#' Example of reference gene expression datasets
#'
#' refData contains processed gene expression data from four datasets,
#' including Systemic Lupus Erythematosus patients and healthy controls.
#' Raw data for each dataset were downloaded from NCBI GEO
#' (GSE65391, GSE45291, GSE61635, and GSE72509, respectively).
#' Platform-dependent preprocessing was performed following established
#' guidelines (Martorell-Marugán et al., 2021). Gene expression data were
#' log2-transformed, and probe sets were annotated to gene symbols.
#' To reduce computational cost in examples, 20 patient and 10 control
#' samples were randomly selected from each dataset.
#'
#' @docType data
#'
#' @usage data(refData)
#'
#' @format An object of class \code{"list"} containing eight objects
#' (dataset1-4 and metadata1-4). Each dataset is a matrix of normalized gene
#' expression values (genes in rows, samples in columns). Each metadata is a
#' dataframe with two columns: samples and group.
"refData"
