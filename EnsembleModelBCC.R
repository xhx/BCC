#' EnsembleModel
#'
#' Ensemble Learning combining models
#' 
#` Author: Xiaohui Xie (xhx@ics.uci.edu)
#' @export


require(rms)
require(MASS)
require(survival)
require(predictiveModeling)


setRefClass(Class = "PredictiveModel")
EnsembleModel <- setRefClass(Class  = "EnsembleModel",
                             contains = "PredictiveModel",
                             fields   = c("models","weights","predictions"),
                             methods  = list(
                               
                               initialize = function(...) {
                                 .self$weights <- c(1,0.2,0.2,0.2,0.2,0.2,0.3 )
                                 return(.self)                                 
                               },
                               
                               
                               customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {        
                                 
                                 cat("train model 1 ..."); flush.console()
                                 ens.mini <- EnsembleModelMinFeatures$new()
                                 ens.mini$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 cat("train model 2 ..."); flush.console()
                                 ens.clnc <- EnsembleModelClncFeatures$new()
                                 ens.clnc$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 cat("train model 3 ..."); flush.console()
                                 ens.meta <- EnsembleModelMetageneFeatures$new()
                                 ens.meta$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 cat("train model 4 ..."); flush.console()
                                 ens.mgnn <- EnsembleModelMetageneNN$new()
                                 ens.mgnn$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 cat("train model 5 ..."); flush.console()
                                 ens.secd <- EnsembleModelSecondaryFeatures$new()
                                 ens.secd$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 cat("train model 6 ..."); flush.console()
                                 ens.gene <- EnsembleModelGeneFeatures$new()
                                 ens.gene$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 cat("train model 7 ..."); flush.console()
                                 ens.clgn <- EnsembleModelClncGeneFeatures$new()
                                 ens.clgn$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                 cat("done!\n"); flush.console()
                                 
                                 .self$models <- c(ens.mini=ens.mini,
                                                   ens.clnc=ens.clnc,
                                                   ens.meta=ens.meta,
                                                   ens.mgnn=ens.mgnn,
                                                   ens.secd=ens.secd,
                                                   ens.gene=ens.gene,
                                                   ens.clgn=ens.clgn
                                 )
                                 
                               },
                               
                               
                               customPredict = function(exprData, copyData, clinicalFeaturesData){                                            
                                 pred <- matrix(0, nrow=nrow(clinicalFeaturesData), ncol=length(.self$models))                                           
                                 for (i in 1:length(.self$models)) {
                                   pred[,i] <- .self$models[[i]]$customPredict(exprData, copyData, clinicalFeaturesData)
                                 }                                            
                                 .self$predictions <- pred                                 
                                 pred <- scale(pred)
                                 
                                 # reweight
                                 p <- 0
                                 for (i in 1:ncol(pred)) 
                                   p <- p + pred[,i]*.self$weights[i]
                                 
                                 return(p)  
                               }
                               
                             )  # end of method list
)  # end of class def
