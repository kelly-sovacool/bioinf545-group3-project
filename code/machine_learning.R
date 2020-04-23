library(tidyverse)

#' Train a random forest model on a training dataset
#'
#' @param training_data  data frame with DiseaseClass column as the outcome
#'                       and all other columns as features
#' @param ncores         number of cores to use for parallel processing
#' @param tune_length    number of attempts for tuning parameters
#'
#' @return model
run_rf <- function(training_data,
                   ncores = 8,
                   tune_length = 50) {
    train_control <- caret::trainControl(
        method = "repeatedcv",
        repeats = 5,
        number = 5,
        returnResamp = "final",
        classProbs = TRUE,
        summaryFunction = caret::twoClassSummary,
        indexFinal = NULL,
        savePredictions = TRUE
    )
    pcluster <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(pcluster)
    model <- caret::train(
        DiseaseClass ~ .,
        data = training_data,
        trControl = train_control,
        method = "rf",
        metric = "ROC",
        tuneLength = tune_length
    )
    parallel::stopCluster(pcluster)
    return(model)
}

#' Split a dataset, train a random forest model, and calculate performance
#'
#' @param data data frame with DiseaseClass column as the outcome and
#'             all other columns as features
#' @param model_filename filename to save the best model to as a TSV file
#' @param conf_mat_filename filename to save the confusion matrix to as an Rds object
#' @param seed random seed
#' @param partition fraction of data to use for training
#' @param ncores number of cores to use for parallel processing
#' @param tune_length number of attempts for tuning parameters
#'
#' @return NULL
predict_rf <- function(data,
                       model_filename,
                       conf_mat_filename,
                       feat_imp_filename,
                       seed = 545,
                       partition = 0.65,
                       ncores = 8,
                       tune_length = 50) {
    set.seed(seed)
    inTraining <- as.vector(caret::createDataPartition(data$DiseaseClass,
                                                       p = partition,
                                                       list = FALSE)
                            )
    training_data <- data[inTraining,]
    testing_data  <- data[-inTraining,]
    rf_model <- run_rf(training_data, ncores, tune_length)
    feat_imp <- caret::varImp(rf_model$finalModel)
    saveRDS(feat_imp, file = feat_imp_filename)
    best_mtry <- rf_model$bestTune$mtry
    best_model <- rf_model$pred %>% filter(mtry == best_mtry)
    readr::write_tsv(best_model, path = model_filename)
    testing_data <- testing_data %>%
        mutate(prediction = predict(rf_model, newdata= testing_data),
               DiseaseClass = as.factor(DiseaseClass))
    conf_mat <- caret::confusionMatrix(data = testing_data$prediction,
                                       reference = testing_data$DiseaseClass,
                                       mode = "everything")
    saveRDS(conf_mat, file = conf_mat_filename)
    return(rf_model)
}
