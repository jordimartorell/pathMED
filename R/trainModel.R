#' Train ML models and perform internal validation
#'
#' @param inputData Numerical matrix or data frame  with samples in columns and
#' features in rows. An ExpressionSet or SummarizedExperiment may also be used.
#' @param metadata Data frame with information for each sample. Samples in rows
#' and variables in columns. If @inputData is an ExpressionSet or
#' SummarizedExperiment, the metadata will be extracted from it.
#' @param models Named list with the ML models generated with
#' caret::caretModelSpec function. methodsML function may be used to prepare
#' this list.
#' @param var2predict Character with the column name of the @metadata to predict
#' @param positiveClass Value that must be considered as positive
#' class (only for categoric variables). If NULL, the last class by
#' alphabetical order is considered as the positive class.
#' @param pairingColumn Optional. Character with the column name of the
#'  @metadata with pairing information (e.g. technical replicates). Paired
#'  samples will always be assigned to the same set (training/test) to avoid
#'  data leakage.
#' @param Koutter Number of outter cross-validation folds.
#' A list of integer with elements for each resampling iteration is admitted.
#' Each list element is a vector of integers corresponding to the rows used for
#' training on that iteration.
#' @param Kinner Number of innter cross-validation folds (for parameter tuning).
#' @param repeatsCV Number of repetitions of the parameter tuning process.
#' @param filterFeatures "rfe" (Recursive Feature Elimination), "sbf" (Selection
#' By Filtering) or NULL (no feature selection).
#' @param filterSizes Only for filterFeatures = "rfe". A numeric vector of
#' integers corresponding to the number of features that should be retained.
#' @param rerank Only for filterFeatures = "rfe". A boolean indicating if the
#' variable importance must be re-calculated each time features are removed.
#' @param continue_on_fail Whether or not to continue training the models if any
#'  of them fail.
#' @param saveLogFile Path to a .txt file in which to save error and warning
#' messages.
#'
#' @return A list with four elements. The first one is the model. The second one
#' is a table with different metrics obtained. The third one is a list with the
#' best parameters selected in tuning process. The last element contains data
#' for AUC plots
#'
#' @author Jordi Martorell-Marugán, \email{jordi.martorell@@genyo.es}
#' @author Daniel Toro-Dominguez, \email{danieltorodominguez@@gmail.com}
#'
#' @importFrom magrittr '%>%'
#' @importFrom methods 'is'
#' @importFrom magrittr '%>%'
#' @import stats
#' @import utils
#'
#' @references Toro-Domínguez, D. et al (2022). \emph{Scoring personalized
#' molecular portraits identify Systemic Lupus Erythematosus subtypes and
#' predict individualized drug responses, symptomatology and
#' disease progression}
#'  . Briefings in Bioinformatics. 23(5)
#'
#' @examples
#' data(exampleData, exampleMetadata)
#'
#' scoresExample <- getScores(exampleData, geneSets = "tmod", method = "GSVA")
#'
#' modelsList <- methodsML("svmLinear", outcomeClass = "character")
#'
#' set.seed(123)
#' trainedModel <- trainModel(
#'     inputData = scoresExample,
#'     metadata = exampleMetadata,
#'     var2predict = "Response",
#'     models = modelsList,
#'     Koutter = 2,
#'     Kinner = 2,
#'     repeatsCV = 1
#' )
#'
#' @export
trainModel <- function(inputData,
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
    saveLogFile = NULL) {
    # 1 Checking
    if (is.null(metadata) & !is(inputData, "ExpressionSet") &
        !is(inputData, "SummarizedExperiment")) {
        stop("If inputData is not an ExpressionSet or SummarizedExperiment,
                metadata must be provided.")
    }

    if (is(inputData, "data.frame")) {
        inputData <- as.matrix(inputData)
    }

    if (is(inputData, "ExpressionSet")) {
        inputData <- Biobase::exprs(inputData)
        metadata <- Biobase::pData(inputData)
    }

    if (is(inputData, "SummarizedExperiment")) {
        inputData <- as.matrix(SummarizedExperiment::assay(inputData))
        metadata <- as.data.frame(SummarizedExperiment::colData(inputData))
    }

    filterSizes <- .trainModelChecking(
        var2predict, metadata, filterFeatures,
        filterSizes, saveLogFile, inputData
    )
    # 2 Filter
    resFilter <- .trainModelFilterRaw(
        inputData, metadata, var2predict,
        positiveClass
    )
    inputData <- resFilter$inputData
    metadata <- resFilter$metadata
    positiveClass <- resFilter$positiveClass
    # # 3 Outcome
    resOutcome <- .trainModelOutcomeClass(
        inputData, metadata, var2predict,
        Koutter, Kinner
    )
    outcomeClass <- resOutcome$outcomeClass
    prior <- resOutcome$prior
    # # 4 Filtering models
    models <- .trainModelMethodsFiltering(inputData, outcomeClass, models)
    # # 5 Koutter
    resKoutter <- .trainModelKoutter(
        inputData, metadata, outcomeClass, Koutter,
        pairingColumn
    )
    sampleSets <- resKoutter$sampleSets
    ntest <- resKoutter$ntest
    # # 6 Training models
    resultNested <- .trainModelResultNested(
        sampleSets, inputData, models,
        pairingColumn, metadata, Kinner,
        repeatsCV, outcomeClass,
        positiveClass,
        filterFeatures, continue_on_fail,
        saveLogFile, rerank, filterSizes
    )
    # 7 Final Features
    Finalfeatures <- .trainModelFinalFeatures(resultNested)
    # 8 Valid Models
    resModels <- .trainModelValidModels(resultNested)
    resultNested <- resModels$resultNested
    validModels <- resModels$validModels
    nullmodels <- resModels$nullmodels
    # 9 perc koutter
    resKoutter <- .trainModelPercKoutter(
        nullmodels, resultNested, validModels,
        models
    )
    validModels <- resKoutter$validModels
    models <- resKoutter$models
    # 10 metrics
    resMetrics <- .trainModelSelectMetrics(
        inputData, outcomeClass,
        positiveClass
    )
    metrics <- resMetrics$metrics
    type <- resMetrics$type
    levels <- resMetrics$levels
    # 11 stats
    stats <- .trainModelStats(
        models, sampleSets, levels, type, metrics,
        resultNested, continue_on_fail, saveLogFile,
        outcomeClass
    )
    # 12 modify stats
    stats <- .trainModelModifyStats(stats, models, prior)
    # 13 parameters, prediction and loss
    resStats <- .trainModelFinalStats(stats, resultNested, ntest)
    stats <- resStats$stats
    parameters <- resStats$parameters
    predsTable <- resStats$predsTable
    # 14 best tune
    bestTune <- .trainModelBestTune(parameters)
    # 15 final model
    fit.model <- .trainModelFit(
        inputData, Finalfeatures, stats, outcomeClass,
        bestTune
    )
    return(list(
        model = fit.model, stats = stats, bestTune = bestTune,
        subsample.preds = predsTable
    ))
}

.trainModelChecking <- function(var2predict, metadata, filterFeatures,
    filterSizes, saveLogFile, inputData) {
    if (!is.null(saveLogFile)) {
        if (file.exists(saveLogFile)) {
            file.remove(saveLogFile) # reset log file
        }
    }

    if (!var2predict %in% colnames(metadata)) {
        stop("var2predict must be a column of metadata")
    }

    if (!is.null(filterFeatures)) {
        if (!filterFeatures %in% c("sbf", "rfe")) {
            stop("filterFeatures must be 'sbf', 'rfe', or NULL")
        } else {
            if (filterFeatures == "rfe") {
                filterSizes <- filterSizes[filterSizes <= nrow(inputData)]
            }
        }
    }
    return(filterSizes)
}

.trainModelFilterRaw <- function(inputData, metadata, var2predict,
    positiveClass) {
    samples <- intersect(rownames(metadata), colnames(inputData))
    if (length(samples) < 1) {
        stop("Row names of metadata and column names of inputData do not match")
    }
    if (is.character(metadata[, var2predict]) |
        is.factor(metadata[, var2predict])) {
        metadata[, var2predict] <- stringi::stri_replace_all_regex(
            metadata[, var2predict],
            pattern = c("/", " ", "-"), replacement = c(".", ".", "."),
            vectorize = FALSE
        )
    }
    if (!is.null(positiveClass)) {
        positiveClass <- stringi::stri_replace_all_regex(
            positiveClass,
            pattern = c("/", " ", "-"),
            replacement = c(".", ".", "."), vectorize = FALSE
        )
    } else {
        positiveClass <- sort(unique(metadata[, var2predict]),
            decreasing = TRUE
        )[1]
    }
    inputData <- inputData[, samples]
    metadata <- metadata[samples, , drop = FALSE]
    # Remove features with all 0
    inputData <- inputData[rowSums(inputData) != 0, ]
    inputData <- data.frame(
        "group" = metadata[, var2predict],
        as.data.frame(t(inputData))
    )
    colnames(inputData) <- stringi::stri_replace_all_regex(
        colnames(inputData),
        pattern = c("/", " ", "-", ":"),
        replacement = c(".", ".", ".", "."), vectorize = FALSE
    )
    inputData <- inputData[!is.na(inputData$group), ]
    if (is.factor(inputData$group)) {
        inputData$group <- as.character(inputData$group)
    }
    if (inputData %>% dplyr::summarise(dplyr::across(
        dplyr::everything(),
        ~ any(is.na(.) |
            is.infinite(.))
    )) %>% any()) {
        stop("There are NA and/or Infinite values in your data, please replace
        or remove them before running trainModel")
    }
    return(list(
        "inputData" = inputData,
        "metadata" = metadata,
        "positiveClass" = positiveClass
    ))
}

.trainModelOutcomeClass <- function(inputData, metadata, var2predict, Koutter,
    Kinner) {
    outcomeClass <- class(inputData$group)
    if (methods::is(inputData$group, "character")) {
        prior <- "MCC"
    } else {
        prior <- "Corr"
    }

    if (outcomeClass == "character" & !is.list(Koutter)) {
        size.m <- min(table(metadata[, var2predict]))
        Kint.m <- as.integer(size.m - (size.m / Koutter))
        if (Kint.m < 2) {
            Kint.m <- 2
        }
        if (size.m < 3) {
            stop("The smallest group has too few samples,
            minimum number of samples per group is 3.")
        }
        if (size.m < 10) {
            message("The smallest group has too few samples,
                the models may not work properly.")
        }
        if (Koutter < 2) {
            stop("Koutter must be 2 or more")
        }
        if (Kinner < 2) {
            stop("Kinner must be 2 or more")
        }
        if (Koutter > size.m) {
            stop("Koutter must be smaller than or equal to the smallest
            group. Maximum Koutter is ", size.m)
        }
        if (Kinner > Kint.m) {
            stop(
                "Not enough samples for ", Koutter,
                " Koutter and ", Kinner, " Kinner. For ", Koutter,
                " Koutter, Kinner must be less or equal to ",
                Kint.m
            )
        }
    }
    return(list(
        "outcomeClass" = outcomeClass,
        "prior" = prior
    ))
}

.trainModelMethodsFiltering <- function(inputData, outcomeClass, models) {
    if (length(unique(inputData$group)) > 2 &
        methods::is(outcomeClass, "character") &
        any(names(models) %in% c("glm", "ada", "gamboost"))) {
        message("glm, ada and gamboost models are not available for
            multi-class models")
        models <- models[!names(models) %in% c("glm", "ada", "gamboost")]
    }
    if (length(models) == 0) {
        stop("No algorithm suitable. Please, check the methodsML function")
    }
    return(models)
}


.trainModelKoutter <- function(inputData, metadata, outcomeClass, Koutter,
    pairingColumn) {
    if (!is.list(Koutter)) {
        if (!is.null(pairingColumn)) {
            isPaired <- metadata[rownames(inputData), pairingColumn]
        } else {
            isPaired <- NULL
        }
        sampleSets <- .makeClassBalancedFolds(
            y = inputData$group, kfold = Koutter,
            repeats = 1, varType = outcomeClass,
            paired = isPaired
        )
    } else {
        sampleSets <- Koutter
    }

    ntest <- sum(unlist(lapply(sampleSets, function(n) {
        nrow(inputData[-as.numeric(unlist(n)), ])
    })))
    return(list(
        "sampleSets" = sampleSets,
        "ntest" = ntest
    ))
}

.trainModelFinalFeatures <- function(resultNested) {
    Finalfeatures <- unique(c(unlist(lapply(resultNested, function(f.ns) {
        tmp.ff <- colnames(f.ns$models[[1]]$trainingData)
        return(tmp.ff[!tmp.ff %in% ".outcome"])
    }))))
}

.trainModelEvaluateModel <- function(resultNested) {
    for (s in seq_len(length(resultNested))) {
        m_toremove <- c()
        for (m in seq_len(length(resultNested[[s]][["models"]]))) {
            if (is.null(resultNested[[s]][["models"]][[m]]) |
                is.null(resultNested[[s]][["preds"]][[m]])) {
                m_toremove <- append(m_toremove, m)
            }
        }
        if (!is.null(m_toremove)) {
            resultNested[[s]][["models"]][[m_toremove]] <- NULL
            resultNested[[s]][["preds"]][[m_toremove]] <- NULL
            resultNested[[s]][["cm"]][[m_toremove]] <- NULL
        }
    }
    return(resultNested)
}

.trainModelNull <- function(resultNested, validModels) {
    nullmodels <- do.call("cbind", lapply(
        seq_len(length(resultNested)),
        function(x) {
            nullmodel <- data.frame(
                row.names = validModels,
                !validModels %in%
                    names(resultNested[[x]][["models"]])
            )
            nullmodel[!nullmodel[, 1], ] <- 0
            colnames(nullmodel) <- paste0("sampleset", x)
            return(nullmodel)
        }
    ))
    return(nullmodels)
}

.trainModelValidModels <- function(resultNested) {
    Finalfeatures <- .trainModelFinalFeatures(resultNested)
    resultNested <- .trainModelEvaluateModel(resultNested)
    validModels <- list()
    validModels <- lapply(resultNested, function(m) {
        vm <- append(validModels, names(m$models))
        return(vm)
    })
    validModels <- unique(unlist(validModels))
    nullmodels <- .trainModelNull(resultNested, validModels)
    return(list(
        "resultNested" = resultNested,
        "nullmodels" = nullmodels,
        "validModels" = validModels
    ))
}


.trainModelResultNested <- function(sampleSets, inputData, models,
    pairingColumn, metadata, Kinner, repeatsCV,
    outcomeClass, positiveClass, filterFeatures,
    continue_on_fail, saveLogFile, rerank, filterSizes) {
    message("Training models...")
    pb <- txtProgressBar(
        min = 0, max = length(sampleSets) * length(models),
        style = 3
    )
    resultNested <- lapply(sampleSets, function(x) {
        training <- inputData[as.numeric(unlist(x)), ]
        testing <- inputData[-as.numeric(unlist(x)), ]

        if (!is.null(pairingColumn)) {
            isPaired <- metadata[rownames(training), pairingColumn]
        } else {
            isPaired <- NULL
        }

        folds <- .makeClassBalancedFolds(
            y = training$group, kfold = Kinner,
            repeats = repeatsCV,
            varType = outcomeClass,
            paired = isPaired
        )

        if (!is.null(filterFeatures)) {
            if (filterFeatures == "sbf") {
                filterCtrl <- caret::sbfControl(
                    functions = NULL, method = "cv",
                    verbose = FALSE,
                    returnResamp = "final",
                    index = folds, allowParallel = TRUE
                )
            }
            if (filterFeatures == "rfe") {
                filterCtrl <- caret::rfeControl(
                    functions = NULL, method = "cv",
                    verbose = FALSE,
                    returnResamp = "final",
                    index = folds, allowParallel = TRUE,
                    rerank = rerank
                )
            }

            tmp <- as.data.frame(training[, 2:ncol(training)])

            if (outcomeClass == "character") {
                y <- factor(training$group)
                if (filterFeatures == "sbf") {
                    filters <- caret::sbf(
                        x = tmp, y = factor(training$group),
                        sbfControl = filterCtrl
                    )
                }
                if (filterFeatures == "rfe") {
                    filters <- caret::rfe(
                        x = tmp, y = factor(training$group, ),
                        rfeControl = filterCtrl,
                        sizes = filterSizes
                    )
                }
            } else {
                y <- as.numeric(training$group)
                if (filterFeatures == "sbf") {
                    filters <- caret::sbf(
                        x = tmp, y = training$group,
                        sbfControl = filterCtrl
                    )
                }
                if (filterFeatures == "rfe") {
                    filters <- caret::rfe(
                        x = tmp, y = training$group,
                        rfeControl = filterCtrl,
                        sizes = filterSizes
                    )
                }
            }


            if (length(filters$optVariables) > 0) {
                training <- training[, c("group", filters$optVariables)]
                testing <- testing[, c("group", filters$optVariables)]
            }
        }

        my_control <- caret::trainControl(
            method = "cv", number = Kinner,
            savePredictions = "final",
            classProbs = ifelse(
                outcomeClass == "character",
                TRUE, FALSE
            ),
            index = folds,
            search = "random"
        )
        global_args <- list(group ~ ., training)
        global_args[["trControl"]] <- my_control

        modelList <- lapply(models, function(m) {
            model_args <- c(global_args, m)
            if (continue_on_fail == TRUE) {
                warn <- err <- NULL
                model <- .removeOutText(
                    withCallingHandlers(
                        tryCatch(do.call(caret::train, model_args),
                            error = function(e) {
                                err <- conditionMessage(e)
                                NULL
                            },
                            warning <- function(w) {
                                warn <<- append(warn, conditionMessage(w))
                                invokeRestart("muffleWarning")
                            }
                        )
                    )
                )

                if (!is.null(saveLogFile) & (!is.null(warn) | !is.null(err))) {
                    write.table(
                        paste0(
                            "Model ", m$method, ":\n",
                            paste0(err, collapse = "\n"),
                            paste0(warn, collapse = "\n"), "\n\n"
                        ),
                        file = saveLogFile,
                        append = TRUE,
                        col.names = FALSE, row.names = FALSE, quote = FALSE
                    )
                }
            } else {
                model <- .removeOutText(do.call(caret::train, model_args))
            }
            setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
            return(model)
        })

        names(modelList) <- names(models)
        nulls <- vapply(modelList, is.null, FUN.VALUE = logical(1))
        modelList <- modelList[!nulls]

        if (length(modelList) == 0) {
            stop("caret:train failed for all models. Please check your data.")
        }

        class(modelList) <- c("caretList")
        modelResults <- .removeOutText(modelList)

        ## Get model stats for Koutter
        predictionTable <- list()
        cm <- list()
        for (m in seq_len(length(modelResults))) {
            if (methods::is(outcomeClass, "character")) {
                classLabels <- levels(as.factor(training$group))
                predTest <- stats::predict(modelResults[[m]],
                    newdata = testing,
                    type = "prob"
                )[, classLabels]
                sel <- unlist(lapply(seq_len(nrow(predTest)), function(n) {
                    ifelse(sum(!is.na(predTest[n, ])) == 0, FALSE,
                        TRUE
                    )
                }))
                if (!any(sel)) {
                    predTest <- NULL
                    cm <- append(cm, list(NULL))
                    names(cm)[m] <- names(modelResults)[m]
                } else {
                    y <- factor(testing$group[sel],
                        levels = unique(testing$group)
                    )
                    predTest <- predTest[sel, ]
                    x <- factor(unlist(lapply(
                        seq_len(nrow(predTest)),
                        function(n) {
                            names(which.max(
                                predTest[n, ]
                            ))
                        }
                    )), levels = unique(testing$group))
                    cmModel <- caret::confusionMatrix(x, y,
                        positive = positiveClass
                    )
                    cm <- append(cm, list(cmModel))
                    names(cm)[m] <- names(modelResults)[m]
                    predTest <- data.frame(predTest, obs = y)
                    rownames(predTest) <- rownames(testing)
                }
            } else { ## Continuous variable
                predTest <- stats::predict(modelResults[[m]], newdata = testing)
                predTest <- data.frame(
                    pred = as.numeric(predTest),
                    obs = as.numeric(testing$group)
                )
                cm <- NULL
            }
            predictionTable <- append(
                predictionTable,
                list(stats::na.omit(predTest))
            )
            names(predictionTable)[m] <- names(modelResults)[m]
        }
        gc()
        return(list(models = modelResults, preds = predictionTable, cm = cm))
    })
    close(pb)
    message("Done")
    return(resultNested)
}

.trainModelFinalFeatures <- function(resultNested) {
    Finalfeatures <- unique(c(unlist(lapply(resultNested, function(f.ns) {
        tmp.ff <- colnames(f.ns$models[[1]]$trainingData)
        return(tmp.ff[!tmp.ff %in% ".outcome"])
    }))))
}

.trainModelFilterModels <- function(resultNested) {
    for (s in seq_len(length(resultNested))) {
        m_toremove <- c()
        for (m in seq_len(length(resultNested[[s]][["models"]]))) {
            if (is.null(resultNested[[s]][["models"]][[m]]) |
                is.null(resultNested[[s]][["preds"]][[m]])) {
                m_toremove <- append(m_toremove, m)
            }
        }
        if (!is.null(m_toremove)) {
            resultNested[[s]][["models"]][[m_toremove]] <- NULL
            resultNested[[s]][["preds"]][[m_toremove]] <- NULL
            resultNested[[s]][["cm"]][[m_toremove]] <- NULL
        }
    }

    validModels <- list()
    validModels <- lapply(resultNested, function(m) {
        vm <- append(validModels, names(m$models))
        return(vm)
    })
    validModels <- unique(unlist(validModels))
    nullmodels <- do.call("cbind", lapply(
        seq_len(length(resultNested)),
        function(x) {
            nullmodel <- data.frame(
                row.names = validModels,
                !validModels %in%
                    names(resultNested[[x]][["models"]])
            )
            nullmodel[!nullmodel[, 1], ] <- 0
            colnames(nullmodel) <- paste0("sampleset", x)
            return(nullmodel)
        }
    ))
    return(list(
        "resultNested" = resultNested,
        "validModels" = validModels,
        "nullmodels" = nullmodels
    ))
}

.trainModelPercKoutter <- function(nullmodels, resultNested, validModels,
                                    models, perc_Koutter_nullModel = 50) {
    removeModel <- rowSums(nullmodels) > length(resultNested) *
        (1 - (perc_Koutter_nullModel / 100))
    removeModel <- names(removeModel[removeModel])
    for (i in seq_len(length(resultNested))) {
        resultNested[[i]][["models"]][names(resultNested[[i]][["models"]])
        %in% removeModel] <- NULL
        resultNested[[i]][["preds"]][names(resultNested[[i]][["preds"]])
        %in% removeModel] <- NULL
    }
    validModels <- validModels[!validModels %in% removeModel]
    failedModels <- names(models)[!names(models) %in% validModels]
    if (length(failedModels) > 0) {
        message(paste0("The following models failed: ", paste0(failedModels,
            collapse = ", "
        )))
    }
    if (is.null(validModels)) {
        stop("All models failed. Please check your data.")
    }
    models <- models[validModels]
    return(list(
        "models" = models,
        "validModels" = validModels
    ))
}

.trainModelSelectMetrics <- function(inputData, outcomeClass, positiveClass) {
    if (outcomeClass == "character") {
        metrics <- c(
            "mcc", "balacc", "accuracy", "recall",
            "specificity", "npv", "precision", "fscore"
        )
        type <- "classification"
        levels <- c(positiveClass, unique(unique(inputData$group))[
            !unique(inputData$group) %in% positiveClass
        ])
    } else {
        metrics <- c("r", "RMSE", "R2", "MAE", "RMAE", "RSE")
        type <- "regression"
        levels <- NULL
    }
    return(list(
        "metrics" = metrics,
        "type" = type,
        "levels" = levels
    ))
}

.trainModelStats <- function(models, sampleSets, levels, type, metrics,
    resultNested, continue_on_fail, saveLogFile,
    outcomeClass) {
    message("Calculating performance metrics...")
    stats <- do.call("cbind", lapply(names(models), function(x) {
        sum.model.x <- do.call("cbind", lapply(
            seq_len(length(sampleSets)),
            function(it) {
                tmp <- resultNested[[it]]$preds[[x]]
                if (methods::is(outcomeClass, "character") & !is.null(tmp)) {
                    obs <- factor(tmp$obs, levels = levels)
                    lab.pred <- unlist(lapply(
                        seq_len(nrow(tmp)),
                        function(n) {
                            names(which.max(tmp[n, !colnames(tmp) %in%
                                "obs"]))
                        }
                    ))
                    lab.pred <- factor(lab.pred, levels = levels)
                } else {
                    obs <- tmp$obs
                    lab.pred <- tmp$pred
                }


                if (continue_on_fail == TRUE) {
                    warn <- err <- NULL
                    resultsTable <- withCallingHandlers(
                        tryCatch(
                            metrica::metrics_summary(
                                obs = obs,
                                pred = lab.pred, type = type,
                                pos_level = 1,
                                metrics_list = metrics
                            ),
                            error = function(e) {
                                err <- conditionMessage(e)
                                NULL
                            }
                        ),
                        warning = function(w) {
                            warn <<- append(warn, conditionMessage(w))
                            invokeRestart("muffleWarning")
                        }
                    )
                    if (!is.null(saveLogFile) &
                        (!is.null(warn) | !is.null(err))) {
                        write.table(
                            paste0(
                                "Model ", x, ", metrics_summary:\n",
                                paste0(err, collapse = "\n"),
                                paste0(warn, collapse = "\n"), "\n\n"
                            ),
                            file = saveLogFile, append = TRUE,
                            col.names = FALSE, row.names = FALSE, quote = FALSE
                        )
                    }
                } else {
                    resultsTable <- suppressWarnings(
                        metrica::metrics_summary(
                            obs = obs,
                            pred = lab.pred,
                            type = type,
                            pos_level = 1,
                            metrics_list = metrics
                        )
                    )
                }
                if (is.null(resultsTable)) {
                    res <- data.frame(Score = rep(NA, length(metrics)))
                } else {
                    rownames(resultsTable) <- resultsTable$Metric
                    res <- data.frame(Score = resultsTable[metrics, "Score"])
                }
                rownames(res) <- metrics
                return(res)
            }
        ))
        sum.model.all <- data.frame(results = apply(
            sum.model.x, 1,
            function(metr) {
                mean(metr, na.rm = TRUE)
            }
        ))
        return(sum.model.all)
    }))

    names(stats) <- names(models)
    return(stats)
}

.trainModelModifyStats <- function(stats, models, prior) {
    if (sum(is.na(stats)) > 0) {
        message("Some metrics could not be calculated and are returned as 0.0")
        stats <- apply(stats, 2, function(x) {
            if (sum(is.na(x)) > 0) {
                corrCol <- x
                corrCol[is.na(corrCol)] <- 0.0
                return(corrCol)
            } else {
                return(x)
            }
        })
    }
    if (length(models) > 1) {
        switch(prior,
            "MCC" = {
                stats <- stats[, order(as.numeric(stats["mcc", ]),
                    decreasing = TRUE
                )]
            },
            "Corr" = {
                stats <- stats[, order(as.numeric(stats["r", ]),
                    decreasing = TRUE
                )]
            }
        )
    }
    return(stats)
}

.trainModelFinalStats <- function(stats, resultNested, ntest) {
    parameters <- do.call("rbind", lapply(resultNested, function(x) {
        x$models[[colnames(stats)[1]]]$bestTune
    }))
    predsTable <- do.call("rbind", lapply(resultNested, function(x) {
        x$preds[colnames(stats)[1]][[1]]
    }))
    lossSamples <- c(unlist(lapply(seq_len(ncol(stats)), function(n) {
        ls.mod <- sum(unlist(lapply(resultNested, function(m) {
            nrow(m$preds[colnames(stats)[n]][[1]])
        })))
    })))
    stats <- rbind(stats, "perc.lossSamples" = 100 - ((lossSamples / ntest) *
        100))
    message("Done")
    return(list(
        "stats" = stats,
        "parameters" = parameters,
        "predsTable" = predsTable
    ))
}

.trainModelBestTune <- function(parameters) {
    bestTune <- data.frame(lapply(seq_len(ncol(parameters)), function(x) {
        tmpValues <- data.frame(parameters[, x])
        if (!methods::is(tmpValues[, 1], "numeric")) {
            tmpValues <- names(table(tmpValues))[order(table(tmpValues),
                decreasing = TRUE
            )][1]
        } else {
            tmpValues <- apply(tmpValues, 2, median)
        }
    }))
    colnames(bestTune) <- paste0(".", colnames(parameters))
    rownames(bestTune) <- NULL
    if (".max_depth" %in% colnames(bestTune)) {
        bestTune$.max_depth <- round(bestTune$.max_depth)
    }
    return(bestTune)
}

.trainModelFit <- function(inputData, Finalfeatures, stats, outcomeClass,
    bestTune) {
    message("Training final model with all samples...")
    newData <- inputData[, colnames(inputData) %in% c("group", Finalfeatures)]

    ## Add aditional parameters for nnet
    if (colnames(stats)[1] == "nnet" & outcomeClass == "character") {
        maxNW <- ((length(Finalfeatures) + 1) * bestTune$.size[1]) +
            ((bestTune$.size[1] + 1) * length(unique(newData$group)))

        fit.model <- withCallingHandlers(tryCatch(
            .removeOutText(
                caret::train(group ~ .,
                    data = newData,
                    method = colnames(stats)[1],
                    tuneGrid = bestTune,
                    MaxNWts = round(maxNW * 10, digits = 0)
                )
            ),
            error = function(e) {
                paste0(
                    "Error fitting the best model (",
                    colnames(stats)[1],
                    ") in all samples. NULL model returned. Try
                            manually selecting a subset of samples and use the
                            optimal parameters provided."
                )
                NULL
            }
        ))
    } else {
        fit.model <- withCallingHandlers(tryCatch(
            .removeOutText(
                caret::train(group ~ .,
                    data = newData,
                    method = colnames(stats)[1],
                    trControl = caret::trainControl(method = "none"),
                    tuneGrid = bestTune
                )
            ),
            error = function(e) {
                paste0(
                    "Error fitting the best model (",
                    colnames(stats)[1],
                    ") in all samples. NULL model returned. Try
                            manually selecting a subset of samples and use the
                            optimal parameters provided."
                )
                NULL
            }
        ))
    }
    message("Done")
    return(fit.model)
}
