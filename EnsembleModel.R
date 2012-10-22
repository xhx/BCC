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
                             fields   = c("idList","models","weights","predictions"),
                             methods  = list(
                               
                               initialize = function(...) {
                                 return(.self)                                 
                               },
                               
                               loadModels = function(idList) {
                                 .self$idList <- idList
                                 
                                 # load models
                                 .self$models <- NULL
                                 for (i in 1:length(.self$idList)) {
                                   cat(i,"\t",.self$idList[i],"\n")
                                   .self$models <- c(.self$models, loadEntity(.self$idList[i])$objects$trainedModel)  
                                 }
                                 
                                 .self$weights <- rep(1, length(.self$models))  # init weights uniformaly
                               },
                               
                               # learn weight according to pred associated with each model
                               setModelWeights = function(clinicalSurvData, pred) {
                                 pred <- scale(pred)
                                 # learn weights based on cox model
                                 #.self$weights = coxph(clinicalSurvData~.,data=as.data.frame(pred))$coef 
                                 
                                 mci <- CIBoost$new(maxIter=6,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                 mci$gridWeight <- ppoints(50)*0.25
                                 mci$customTrain(pred, clinicalSurvData,initWithMaxCI=TRUE)  #init with max CI genes
                                 w <- rep(0, length(.self$models))
                                 for (i in 1:nrow(mci$weights)) {
                                   w[mci$weights[i,1]] <- w[mci$weights[i,1]] + mci$weights[i,2]
                                 }
                                 .self$weights <- w
                                 .self$wci <- mci
                               },
                               
                               
                               customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {                                           
                                 #pred <- matrix(0, nrow=nrow(clinicalFeaturesData), ncol=length(.self$idList))                                            
                                 
                                 for (i in 1:length(.self$idList)) {
                                   cat("train model ", i,"\t",.self$idList[i],"\n")                                             
                                   .self$models[[i]]$customTrain(exprData, copyData, clinicalFeaturesData, clinicalSurvData)
                                   #pred[,i] <- .self$models[[i]]$customPredict(exprData, copyData, clinicalFeaturesData)
                                 }
                                 
                                 #pred <- scale(pred)
                                 #setModelWeights(clinicalSurvData, pred)                                                                                    
                               },
                               
                               
                               customPredict = function(exprData, copyData, clinicalFeaturesData){                                            
                                 pred <- matrix(0, nrow=nrow(clinicalFeaturesData), ncol=length(.self$idList))                                           
                                 for (i in 1:length(.self$idList)) {
                                   pred[,i] <- .self$models[[i]]$customPredict(exprData, copyData, clinicalFeaturesData)
                                 }                                            
                                 .self$predictions <- pred
                                 
                                 pred <- scale(pred)
                                 
                                 p <- 0
                                 for (i in 1:ncol(pred)) 
                                   p <- p + pred[,i]*.self$weights[i]
                                 
                                 return(p)  
                               }
                               
                             )  # end of method list
)  # end of class def
