% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mScores_filterPaths.R
\name{mScores_filterPaths}
\alias{mScores_filterPaths}
\title{Filter pathways from the reference M-scores dataset}
\usage{
mScores_filterPaths(
  MRef,
  min_datasets = round(length(MRef[[1]]) * 0.34),
  perc_samples = 10,
  Pcutoff = 0.05,
  plotMetrics = TRUE
)
}
\arguments{
\item{MRef}{output from the mScores_createReference function}

\item{min_datasets}{number of datasets that each pathway must meet the
perc_samples threshold}

\item{perc_samples}{minimun percentage of samples in a dataset in which a
pathway must be significant}

\item{Pcutoff}{P-value cutoff for significance}

\item{plotMetrics}{Plot number of significant pathways selected based on the
different combination of perc_samples and min_datasets parameters}
}
\value{
A list with the selected pathways
}
\description{
Filter pathways from the reference M-scores dataset
}
\examples{
data(refData)

refObject <- buildRefObject(
    data = list(
        refData$dataset1, refData$dataset2,
        refData$dataset3, refData$dataset4
    ),
    metadata = list(
        refData$metadata1, refData$metadata2,
        refData$metadata3, refData$metadata4
    ),
    groupVar = "group",
    controlGroup = "Healthy_sample"
)

exampleRefMScore <- mScores_createReference(refObject, geneSets = "tmod")
relevantPaths <- mScores_filterPaths(exampleRefMScore, min_datasets = 3)
}
\references{
Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
molecular portraits identify Systemic Lupus Erythematosus subtypes and
predict individualized drug responses, symptomatology and
disease progression}
 . Briefings in Bioinformatics. 23(5)
}
\seealso{
\code{\link{mScores_createReference}}
}
\author{
Jordi Martorell-Marugán, \email{jordi.martorell@genyo.es}

Daniel Toro-Dominguez, \email{danieltorodominguez@gmail.com}
}
