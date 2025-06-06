% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trainModel.R
\name{trainModel}
\alias{trainModel}
\title{Train ML models and perform internal validation}
\usage{
trainModel(
  inputData,
  metadata = NULL,
  models = methodsML(outcomeClass = "character"),
  var2predict,
  positiveClass = NULL,
  pairingColumn = NULL,
  Koutter = 5,
  Kinner = 4,
  repeatsCV = 5,
  filterFeatures = NULL,
  filterSizes = seq(2, 100, by = 2),
  rerank = FALSE,
  continue_on_fail = TRUE,
  saveLogFile = NULL,
  use.assay = 1
)
}
\arguments{
\item{inputData}{Numerical matrix or data frame  with samples in columns and
features in rows. An ExpressionSet or SummarizedExperiment may also be used.}

\item{metadata}{Data frame with information for each sample. Samples in rows
and variables in columns. If @inputData is an ExpressionSet or
SummarizedExperiment, the metadata will be extracted from it.}

\item{models}{Named list with the ML models generated with
caret::caretModelSpec function. methodsML function may be used to prepare
this list.}

\item{var2predict}{Character with the column name of the @metadata to predict}

\item{positiveClass}{Value that must be considered as positive
class (only for categoric variables). If NULL, the last class by
alphabetical order is considered as the positive class.}

\item{pairingColumn}{Optional. Character with the column name of the
@metadata with pairing information (e.g. technical replicates). Paired
samples will always be assigned to the same set (training/test) to avoid
data leakage.}

\item{Koutter}{Number of outter cross-validation folds.
A list of integer with elements for each resampling iteration is admitted.
Each list element is a vector of integers corresponding to the rows used for
training on that iteration.}

\item{Kinner}{Number of innter cross-validation folds (for parameter tuning).}

\item{repeatsCV}{Number of repetitions of the parameter tuning process.}

\item{filterFeatures}{"rfe" (Recursive Feature Elimination), "sbf" (Selection
By Filtering) or NULL (no feature selection).}

\item{filterSizes}{Only for filterFeatures = "rfe". A numeric vector of
integers corresponding to the number of features that should be retained.}

\item{rerank}{Only for filterFeatures = "rfe". A boolean indicating if the
variable importance must be re-calculated each time features are removed.}

\item{continue_on_fail}{Whether or not to continue training the models if any
of them fail.}

\item{saveLogFile}{Path to a .txt file in which to save error and warning
messages.}

\item{use.assay}{If SummarizedExperiments are used, the number of the assay 
to extract the data.}
}
\value{
A list with four elements. The first one is the model. The second one
is a table with different metrics obtained. The third one is a list with the
best parameters selected in tuning process. The last element contains data
for AUC plots
}
\description{
Train ML models and perform internal validation
}
\examples{
data(pathMEDExampleData, pathMEDExampleMetadata)

scoresExample <- getScores(pathMEDExampleData, geneSets = "tmod", 
                             method = "GSVA")

modelsList <- methodsML("svmLinear", outcomeClass = "character")

set.seed(123)
trainedModel <- trainModel(
    inputData = scoresExample,
    metadata = pathMEDExampleMetadata,
    var2predict = "Response",
    models = modelsList,
    Koutter = 2,
    Kinner = 2,
    repeatsCV = 1
)

}
\references{
Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
molecular portraits identify Systemic Lupus Erythematosus subtypes and
predict individualized drug responses, symptomatology and
disease progression}
 . Briefings in Bioinformatics. 23(5)
}
\author{
Jordi Martorell-Marugán, \email{jordi.martorell@genyo.es}

Daniel Toro-Dominguez, \email{danieltorodominguez@gmail.com}
}
