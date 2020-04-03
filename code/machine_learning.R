run_rf <- function(data) {
    cross_val <- trainControl(method="repeatedcv",
                          repeats = 1,
                          number=5,
                          returnResamp="final",
                          classProbs=TRUE,
                          summaryFunction=twoClassSummary,
                          indexFinal=NULL,
                          savePredictions = TRUE)
    model <- train(DiseaseClass ~ .,
               data = data,
               trControl = cross_val,
               method="rf",
               metric="ROC",
               tuneLength=5)
    return(model)
}
