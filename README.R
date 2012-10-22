# Models for Sage BCC
# Instructions on how to use each model
#

modelClassFile <- "EnsembleModelBCC.R"
modelClassFile2 <- "EnsembleModel.R"
suppFuncFile <- "EnsembleSuppFunc.R"


#1 ENSEMBLE 101501
source(suppFuncFile)
predModel <- source(modelClassFile)$value$new()
predModel$customTrain(trainingData$exprData, trainingData$copyData,trainingData$clinicalFeaturesData, trainingData$clinicalSurvData)
predModel$weights <- c(1,0.2,0,0.2,0.2,0.2,0.3)

#2 ENSEMBLE 1015 SUB
source(suppFuncFile)
predModel <- source(modelClassFile)$value$new()
predModel$customTrain(trainingData$exprData, trainingData$copyData,trainingData$clinicalFeaturesData, trainingData$clinicalSurvData)
predModel$weights <- c(1,0.15,0.15,0.15,0.15,0.15,0.2)

#3 ENSEMBLE 101502
source(suppFuncFile)
predModel <- source(modelClassFile)$value$new()
predModel$customTrain(trainingData$exprData, trainingData$copyData,trainingData$clinicalFeaturesData, trainingData$clinicalSurvData)
predModel$weights <- c(1,0.3,0,0.3,0.3,0.3,0.5)

#4 ENSEMBLE 100902
source(suppFuncFile)
predModel <- source(modelClassFile2)$value$new()
idList <- c("syn1420933", "syn1420549", "syn1420939", "syn1417636")
predModel$loadModels(idList)
predModel$customTrain(trainingData$exprData, trainingData$copyData,trainingData$clinicalFeaturesData, trainingData$clinicalSurvData)
predModel$weights <- c(1,0.2,0.2,0.3)
 
#4 ENSEMBLE 101001
source(suppFuncFile)
predModel <- source(modelClassFile2)$value$new()
idList <- c('syn1435273', 'syn1420549', 'syn1420939', 'syn1417636') 
predModel$loadModels(idList)
predModel$customTrain(trainingData$exprData, trainingData$copyData,trainingData$clinicalFeaturesData, trainingData$clinicalSurvData)
predModel$weights <- c(1,0.2,0.2,0.3)
 
