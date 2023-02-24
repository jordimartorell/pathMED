#' Example of reference gene expression data
#'
#' refData contains processed gene expression data for 4 datasets
#' including lupus and healthy samples. Raw data for each datasets was
#' downloaded from NCBI GEO (GSE65391, GSE45291, GSE61635, GSE72509). Raw data
#' platform-dependent preprocessing was performed following previously
#' guidelines (Martorell-Marugán et al., 2021). Finally, gene expression was
#' log2-transformed and genes were annotated from probe sets to gene symbol. 20
#' patients and 10 control samples were randomly selected with the aim of
#' making the examples less computationally expensive.
#'
#' @docType data
#'
#' @usage data(refData)
#'
#' @format An object of class \code{"list"} with 4 datasets. Each dataset is a
#' list with two data.frames ("Disease" and "Healthy"). Each data frame
#' contains the normalized gene expression, with genes in rows and samples in
#' columns.
#'
"refData"
