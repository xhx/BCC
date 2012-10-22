#'
#' Functions and Classes for BCC challenge. 
#'
#' author: Xiaohui Xie (xhx@ics.uci.edu)
#

require(rms)
require(MASS)
require(survival)
require(predictiveModeling)
require(randomSurvivalForest)
require(rpart)
require(caret)


setRefClass(Class = "PredictiveModel")

if(!require(gbm)){
  download.file("http://cran.r-project.org/src/contrib/gbm_1.6-3.2.tar.gz", destfile="./gbm_1.6-3.2.tar.gz")
  install.packages("./gbm_1.6-3.2.tar.gz", repos=NULL)
  library(gbm)
}

if(!require(DreamBox7)){
  download.file("http://www.ics.uci.edu/~xhx/out/BCC/DreamBox7_0.21.tar.gz", destfile="./DreamBox7_0.21.tar.gz")
  install.packages("./DreamBox7_0.21.tar.gz", repos=NULL)
  library(DreamBox7)
}


ImputeMissingClinicalFeatures <- function(c, random=FALSE) {
  for (i in 1:ncol(c)) {
    if (sum(is.na(c[,i]))==0) next  # no NA in this column
    
    if (is.numeric(c[,i])) {
      mu <- mean(c[,i], na.rm=TRUE)
      sd <- sd(c[,i], na.rm=TRUE)
      if (!random) 
        c[which(is.na(c[,i])),i] <- mu
      else 
        c[which(is.na(c[,i])),i] <- rnorm(sum(is.na(c[,i])),mu,sd)
    } else if (is.factor(c[,i])) {
      foo.levels <- levels(c[,i])
      foo.freq <- tapply(rep(1,length(c[,i])),c[,i],sum)
      foo.freq <- na.omit(as.matrix(foo.freq))
      if (!random)
        c[which(is.na(c[,i])),i] <- rownames(foo.freq)[which.max(foo.freq)]
      else
        c[which(is.na(c[,i])),i] <- sample(rownames(foo.freq),size=sum(is.na(c[,i])),replace=TRUE,prob=foo.freq/sum(foo.freq))
    }
  }
  return(c)
}

GenerateCatClncFeatures <- function(c){
  
  age.cat <- recode(c$age_at_diagnosis, "lo:40=1; 40.0001:75=2; 75.0001:hi=3")
  lym.cat <- recode(c$lymph_nodes_positive, "0=1; 1:5=2; else=3")
  tr.CT <- c$tr.CT
  
  cmat <- data.frame(tr.CT, age.cat, lym.cat)
  return(cmat)
}

# from the Goldi model
# use_osloval: use_non_osloval features?
GoldiExpandClnc <- function(c, use_osloval=TRUE) { 
  h.IDC = as.numeric(c$histological_type == "IDC")
  h.ILC = as.numeric(c$histological_type == "ILC")
  h.IDCpILC = as.numeric(c$histological_type == "IDC+ILC")
  h.IDCnMED = as.numeric(c$histological_type == "IDC-MED")
  h.IDCnMUC = as.numeric(c$histological_type == "IDC-MUC")
  h.IDCnTUB = as.numeric(c$histological_type == "IDC-TUB")
  h.other = as.numeric(c$histological_type == "MIXED NST AND A SPECIAL TYPE" | 
    c$histological_type == "OTHER" | c$histological_type == 
    "OTHER INVASIVE" | c$histological_type == "INVASIVE TUMOR" | 
    c$histological_type == "PHYL")
  er.P = as.numeric(c$ER_IHC_status == "pos")
  er.N = as.numeric(c$ER_IHC_status == "neg")
  tr.CT = as.numeric((c$Treatment == "CT") | (c$Treatment == 
    "CT/HT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == 
    "CT/RT"))
  tr.HT = as.numeric((c$Treatment == "HT") | (c$Treatment == 
    "CT/HT") | (c$Treatment == "HT/RT") | (c$Treatment == 
    "CT/HT/RT"))
  tr.RT = as.numeric((c$Treatment == "RT") | (c$Treatment == 
    "CT/HT/RT") | (c$Treatment == "CT/RT") | (c$Treatment == 
    "HT/RT"))
  trMat = cbind(tr.CT, tr.HT, tr.RT)
  gd.1 = as.numeric(c$grade == 1)
  gd.2 = as.numeric(c$grade == 2)
  gd.3 = as.numeric(c$grade == 3)
  
  her2.snp6.gain = as.numeric(c$HER2_SNP6_state == "GAIN")
  her2.snp6.loss = as.numeric(c$HER2_SNP6_state == "LOSS")
  her2.snp6.neut = as.numeric(c$HER2_SNP6_state == "NEUT")
  
  if (use_osloval) {
    grp.1 = as.numeric(c$NOT_IN_OSLOVAL_group == 1)
    grp.2 = as.numeric(c$NOT_IN_OSLOVAL_group == 2)
    grp.3 = as.numeric(c$NOT_IN_OSLOVAL_group == 3)
    grp.4 = as.numeric(c$NOT_IN_OSLOVAL_group == 4)
    stg.0 = as.numeric(c$NOT_IN_OSLOVAL_stage == 0)
    stg.1 = as.numeric(c$NOT_IN_OSLOVAL_stage == 1)
    stg.2 = as.numeric(c$NOT_IN_OSLOVAL_stage == 2)
    stg.3 = as.numeric(c$NOT_IN_OSLOVAL_stage == 3 | c$NOT_IN_OSLOVAL_stage == 4)  
    st.1 = as.numeric(c$NOT_IN_OSLOVAL_Site == 1)
    st.2 = as.numeric(c$NOT_IN_OSLOVAL_Site == 2)
    st.3 = as.numeric(c$NOT_IN_OSLOVAL_Site == 3)
    st.4 = as.numeric(c$NOT_IN_OSLOVAL_Site == 4)
    st.5 = as.numeric(c$NOT_IN_OSLOVAL_Site == 5)
    
    cmat <- data.frame(c[, c(1:3)], gd.1, gd.2, gd.3, h.IDC, 
                       h.ILC, h.IDCpILC, h.IDCnMED, h.IDCnMUC, h.IDCnTUB, h.other, 
                       er.N, er.P, tr.CT, tr.HT, tr.RT, her2.snp6.gain, her2.snp6.loss, 
                       her2.snp6.neut, grp.1, grp.2, grp.3, grp.4, stg.0, stg.1, 
                       stg.2, stg.3, st.1, st.2, st.3, st.4, st.5)
  } else {
    cmat <- data.frame(c[, c(1:3)], gd.1, gd.2, gd.3, h.IDC, 
                       h.ILC, h.IDCpILC, h.IDCnMED, h.IDCnMUC, h.IDCnTUB, h.other, 
                       er.N, er.P, tr.CT, tr.HT, tr.RT, her2.snp6.gain, her2.snp6.loss, 
                       her2.snp6.neut)
  }
  
  
  for (i in 4:ncol(cmat)) {
    cmat[, i] = factor(cmat[, i])
  }
  return(cmat)
}


GenerateClncFeatures <- function(c){
  
  lym.P <- as.numeric(c$lymph_nodes_positive)
  
  h.IDC =as.numeric(c$histological_type=="IDC")
  h.ILC =as.numeric(c$histological_type=="ILC")
  h.IDCpILC =as.numeric(c$histological_type=="IDC+ILC")
  h.IDCnTUB =as.numeric(c$histological_type=="IDC-TUB")
  h.IDCnMUC =as.numeric(c$histological_type=="IDC-MUC")
  h.IDCnMED =as.numeric(c$histological_type=="IDC-MED")
  h.MIXED =as.numeric(c$histological_type=="MIXED NST AND A SPECIAL TYPE")
  h.OTHER =as.numeric(c$histological_type=="OTHER")
  h.OTHERINV =as.numeric(c$histological_type=="OTHER INVASIVE")
  h.INVTUMOR =as.numeric(c$histological_type=="INVASIVE TUMOR")
  h.others = h.MIXED + h.OTHER + h.OTHERINV + h.INVTUMOR
  
  er.P=as.numeric(c$ER_IHC_status=="pos")
  er.N=as.numeric(c$ER_IHC_status=="neg")
  
  c.l=as.numeric(c$NOT_IN_OSLOVAL_cellularity=="low")
  c.m=as.numeric(c$NOT_IN_OSLOVAL_cellularity=="moderate")
  c.h=as.numeric(c$NOT_IN_OSLOVAL_cellularity=="high")
  
  
  p.LumA=as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype=="LumA")
  p.LumB=as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype=="LumB")                                         
  p.Her2=as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype=="Her2")                                         
  p.Normal=as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype=="Normal")         
  p.Basal=as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype=="Basal") 
  
  tr.CT = as.numeric((c$Treatment == "CT") | (c$Treatment == "CT/HT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == "CT/RT"))
  tr.HT = as.numeric((c$Treatment == "HT") | (c$Treatment == "CT/HT") | (c$Treatment == "HT/RT") | (c$Treatment =="CT/HT/RT"))
  tr.RT = as.numeric((c$Treatment == "RT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == "CT/RT") | (c$Treatment == "HT/RT"))
  
  st.1 = as.numeric(c$NOT_IN_OSLOVAL_Site == 1)
  st.2 = as.numeric(c$NOT_IN_OSLOVAL_Site == 2)
  st.3 = as.numeric(c$NOT_IN_OSLOVAL_Site == 3)
  
  stage <- as.numeric(c$NOT_IN_OSLOVAL_stage)
  age <- as.numeric(c$age_at_diagnosis)
  size <- as.numeric(c$size)
  
  
  
  cmat<-data.frame(stage, age, size, lym.P, h.IDCnTUB, h.IDCnMUC, h.IDCnMED,h.others, h.IDC, h.ILC,h.IDCpILC,
                   er.P, er.N, tr.CT, tr.HT, tr.RT, p.LumA, p.LumB, p.Her2, p.Normal, p.Basal, st.1, st.2, st.3)
  for(i in 5:ncol(cmat)){
    cmat[,i] = factor(cmat[,i])
  }
  return(cmat)
}



bcc.train.gene.cox <- function(exprData, clinicalSurvData, selectedGenes=NULL) {
  # train cox model based on selected genes
  # Args: 
  # Return:
  #
  
  # extract data  
  tmpData <- ExtractGeneExprData(exprData, selectedGenes) 
  
  tmpModel <- coxph(clinicalSurvData ~(.), data = tmpData, model=TRUE)
  #tmpModel <- stepAIC(tmpModel, direction="both")    
  return(tmpModel)
}


bcc.train.gene.gbm <- function(exprData, clinicalSurvData, selectedGenes=NULL) {
  # train gbm model based on selected genes
  # Args: 
  # Return:
  
  # extract data
  tmpData <- ExtractGeneExprData(exprData, selectedGenes) 
  
  tmpModel <- gbm.fit(tmpData,clinicalSurvData,distribution="coxph", shrinkage=0.001, n.trees=1000,
                      interaction.depth=4, bag.fraction=0.8, train.fraction=1, verbose=F)
  return(tmpModel)
}



ExtractGeneExprData <- function(exprData, geneList=NULL) {
  # extract gene expression as data.frame, selected a subset of genes of geneList is specified
  # Args:
  # Return:
  #
  
  if (class(exprData)=="ExpressionSet") {  
    data.expr <- exprs(exprData)
    data.expr[is.na(data.expr)] <- median(data.expr, na.rm=TRUE)
    data.expr <- as.data.frame(t(data.expr))
  } else
    data.expr <- exprData
  
  # if a gene list is specified
  if (is.null(geneList)==FALSE) {
    geneList <- intersect(geneList, colnames(data.expr))
    data.expr  <- data.expr[, geneList]
  }
  return(data.expr)
}



bcc.train.clnc.cox <- function(clinicalFeaturesData,clinicalSurvData, selectedFeatures=NULL) {
  # train cox model based on clinical features
  #
  if(class(clinicalSurvData) != "Surv"){
    stop("Expecting 'responseData' object of type 'Surv'")
  } 
  
  clinical <- clinicalFeaturesData
  
  if (is.null(selectedFeatures)==FALSE)
    clinical <- clinical[,selectedFeatures]
  
  upper = terms(clinicalSurvData~(.), data = clinical)
  coxmodel.clnc = step(coxph(clinicalSurvData~1, data=clinical), scope=upper, direction="both", k=2, trace=FALSE)
  #coxmodel.clnc = coxph(clinicalSurvData~., data=data.frame(clinical))
  
  return(coxmodel.clnc)
}    


bcc.train.clnc.gbm <- function(clinicalFeaturesData,clinicalSurvData,selectedFeatures=NULL) {
  # train gbm model based on  clinical features
  if(class(clinicalSurvData) != "Surv"){
    stop("Expecting 'responseData' object of type 'Surv'")
  } 
  
  clinical <- clinicalFeaturesData
  
  if (is.null(selectedFeatures)==FALSE)
    clinical <- clinical[,selectedFeatures]
  
  gbmmodel.clnc = gbm.fit(clinical,clinicalSurvData,distribution="coxph", shrinkage=0.001, n.trees=2000, interaction.depth=7, bag.fraction=0.8, train.fraction=1, verbose=F)
  
  return(gbmmodel.clnc)
}



bcc.train.expandedClnc.cox <- function(clinicalFeaturesData,clinicalSurvData, selectedFeatures=NULL) {
  # train cox model based on expanded clinical features
  #
  if(class(clinicalSurvData) != "Surv"){
    stop("Expecting 'responseData' object of type 'Surv'")
  } 
  
  clinical <- GenerateClncFeatures(clinicalFeaturesData)
  
  if (is.null(selectedFeatures)==FALSE)
    clinical <- clinical[,selectedFeatures]
  
  upper = terms(clinicalSurvData~(.), data = clinical)
  coxmodel.clnc = step(coxph(clinicalSurvData~1, data=clinical), scope=upper, direction="both", k=2, trace=FALSE)
  #coxmodel.clnc = coxph(clinicalSurvData~., data=data.frame(clinical))
  
  return(coxmodel.clnc)
}    


bcc.train.expandedClnc.gbm <- function(clinicalFeaturesData,clinicalSurvData,selectedFeatures=NULL) {
  # train gbm model based on expanded clinical features
  if(class(clinicalSurvData) != "Surv"){
    stop("Expecting 'responseData' object of type 'Surv'")
  } 
  
  clinical <- GenerateClncFeatures(clinicalFeaturesData)
  
  if (is.null(selectedFeatures)==FALSE)
    clinical <- clinical[,selectedFeatures]
  
  gbmmodel.clnc = gbm.fit(clinical,clinicalSurvData,distribution="coxph", shrinkage=0.001, n.trees=2000, interaction.depth=7, bag.fraction=0.8, train.fraction=1, verbose=F)
  
  return(gbmmodel.clnc)
}



# normalize each gene to be mean 0 and var 1 across different samples
ScaleNormExprData <- function(exprData) {
  # exprData: numSamples x numGenes   matrix
  return(scale(exprData))
} 

# array quantile normalize each gene across different samples
ArrayQuantileNormExprData <- function(exprData) {
  expr <- exprs(exprData)
  # expr: numGenes  x numSamples matrix
  for (j in 1:ncol(expr)) 
    expr[,j] <- MapToNorm(expr[,j])
  
  return(ExpressionSet(expr))
} 



# quantile normalize each gene across different samples
QuantileNormExprData <- function(exprData) {
  # exprData: numSamples x numGenes  matrix
  for (j in 1:ncol(exprData)) 
    exprData[,j] <- MapToNorm(exprData[,j])
  
  return(exprData)
} 

MapToNorm <- function(x) {
  # map data vector x to norm based on quantile 
  # 
  return( qnorm(ppoints(length(x)))[order(order(x))] )
}


# from Warwick systems biology 
ExtractClinicalData <- function(data.cl){
  ##EXTRACT THE BASIC FEATURES
  itemNames         = rownames(data.cl)
  age               = as.numeric(data.cl$age_at_diagnosis)
  young             = as.numeric(age<40)
  old               = as.numeric(age>77)
  size              = as.numeric(data.cl$size)
  lymphNodes        = as.numeric(data.cl$lymph_nodes_positive)
  grade             = as.numeric(data.cl$grade)
  prStatus.expr     = as.numeric(data.cl$PR.Expr=="+")
  her2Status        = as.numeric(data.cl$HER2_IHC_status)-1
  her2.snp6.gain    = as.numeric(data.cl$HER2_SNP6_state=="GAIN")
  her2.snp6.neutral = as.numeric(data.cl$HER2_SNP6_state=="NEUT")
  her2.snp6.loss    = as.numeric(data.cl$HER2_SNP6_state=="LOSS")
  her2Status.expr   = as.numeric(data.cl$Her2.Expr=="+")                         
  ##CONVERT NON-NUMERIC: HISTOLOGICAL TYPE
  ##note that almost 80% of cases are IDC
  ##these features may be too sparse to be useful??
  ##Maybe try histo_idc?
  working              <- data.cl$histological_type
  histo_dcis           <- as.numeric(working=="DCIS")
  histo_phyl           <- as.numeric(working=="PHYL")
  histo_idc            <- as.numeric(working=="IDC")
  histo_idc_ilc        <- as.numeric(working=="IDC+ILC")
  histo_idc_med        <- as.numeric(working=="IDC-MED")
  histo_idc_muc        <- as.numeric(working=="IDC-MUC")
  histo_idc_tub        <- as.numeric(working=="IDC-TUB")
  histo_ilc            <- as.numeric(working=="ILC")
  histo_mixed          <- as.numeric(working=="MIXED NST AND A SPECIAL TYPE")
  histo_other          <- as.numeric(working=="OTHER")
  histo_otherInvasive  <- as.numeric(working=="OTHER INVASIVE")  ##empty feature! (but not in test set)
  histo_invasiveTumour <- as.numeric(working=="INVASIVE TUMOUR") ##empty feature! (but not in test set)
  histo_combo          <- histo_mixed + histo_other + histo_otherInvasive + histo_invasiveTumour
  ##ER STATUS
  ##patch the NAs in this feature
  #erStatus          = as.numeric(data.cl$ER_IHC_status=="neg")
  erStatus.expr     = as.numeric(data.cl$ER.Expr=="+")
  erStatus          = as.numeric(data.cl$ER_IHC_status) - 1
  patch             = is.na(erStatus)
  erStatus[patch]   = erStatus.expr[patch]                         
  ##CONVERT NON-NUMERIC: TREATMENT
  working      <- data.cl$Treatment
  treatment_ct <- as.numeric(working=="CT" | working=="CT/HT" | working=="CT/RT" | working=="CT/HT/RT")
  treatment_ht <- as.numeric(working=="HT" | working=="CT/HT" | working=="HT/RT" | working=="CT/HT/RT")
  treatment_rt <- as.numeric(working=="RT" | working=="HT/RT" | working=="CT/RT" | working=="CT/HT/RT")   
  ##FEATURES THAT ONLY APPLY TO THE METABRIC DATA
  postMenopausal    = as.numeric(data.cl$NOT_IN_OSLOVAL_menopausal_status_inferred=="post")
  stage             = as.numeric(data.cl$NOT_IN_OSLOVAL_stage)
  lymphNodesRemoved = as.numeric(data.cl$NOT_IN_OSLOVAL_lymph_nodes_removed)
  NPI               = as.numeric(data.cl$NOT_IN_OSLOVAL_NPI)
  p53status_wt      = as.numeric(data.cl$NOT_IN_OSLOVAL_P53_mutation_status=="WT")
  p53status_mut     = as.numeric(data.cl$NOT_IN_OSLOVAL_P53_mutation_status=="MUT")
  ##(METABRIC) Genefu
  ##under construction!
  #Genefu            = as.numeric(data.cl$)                         
  ##CURTIS ET AL SUBTYPES
  intClust_1  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="1")
  intClust_2  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="2")
  intClust_3  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="3")
  intClust_4  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="4")
  intClust_5  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="5")
  intClust_6  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="6")
  intClust_7  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="7")
  intClust_8  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="8")
  intClust_9  = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="9")
  intClust_10 = as.numeric(data.cl$NOT_IN_OSLOVAL_IntClustMemb=="10")                         
  ##(METABRIC) P53 INFORMATION
  ##under construction!
  #p53type           = as.numeric(data.cl$)
  #p53details        = as.numeric(data.cl$)                         
  ##(METABRIC) GROUP
  group       = as.numeric(data.cl$NOT_IN_OSLOVAL_group)
  group_1     = as.numeric(data.cl$NOT_IN_OSLOVAL_group=="1")
  group_2     = as.numeric(data.cl$NOT_IN_OSLOVAL_group=="2")
  group_3     = as.numeric(data.cl$NOT_IN_OSLOVAL_group=="3")
  group_4     = as.numeric(data.cl$NOT_IN_OSLOVAL_group=="4")
  group_other = as.numeric(data.cl$NOT_IN_OSLOVAL_group=="other")                         
  ##(METABRIC) PAM 50 SUBTYPE
  working      <- data.cl$NOT_IN_OSLOVAL_Pam50Subtype
  pam50_basal  <- as.numeric(working=="Basal") 
  pam50_her2   <- as.numeric(working=="Her2") 
  pam50_lumA   <- as.numeric(working=="LumA") 
  pam50_lumB   <- as.numeric(working=="LumB") 
  pam50_normal <- as.numeric(working=="Normal") 
  ##(METABRIC) CELLULARITY
  cellularity_low  = as.numeric(data.cl$NOT_IN_OSLOVAL_cellularity=="low")
  cellularity_mod  = as.numeric(data.cl$NOT_IN_OSLOVAL_cellularity=="moderate")
  cellularity_high = as.numeric(data.cl$NOT_IN_OSLOVAL_cellularity=="high")                         
  ##(METABRIC) SITE
  site_1 = as.numeric(data.cl$NOT_IN_OSLOVAL_Site=="1")
  site_2 = as.numeric(data.cl$NOT_IN_OSLOVAL_Site=="2")
  site_3 = as.numeric(data.cl$NOT_IN_OSLOVAL_Site=="3")
  site_4 = as.numeric(data.cl$NOT_IN_OSLOVAL_Site=="4")
  site_5 = as.numeric(data.cl$NOT_IN_OSLOVAL_Site=="5")
  ##CONSTRUCT, RETURN A NEW DATA FRAME
  outputData <- data.frame(age, size, lymphNodes, grade, erStatus,
                           erStatus.expr, prStatus.expr, her2Status,
                           her2.snp6.gain, her2.snp6.neutral, her2.snp6.loss,
                           her2Status.expr, histo_dcis, histo_phyl,
                           histo_idc, histo_idc_ilc, histo_idc_med,
                           histo_idc_muc, histo_idc_tub, histo_ilc,
                           histo_mixed, histo_other, histo_otherInvasive,
                           histo_invasiveTumour, histo_combo,
                           treatment_ct, treatment_ht, treatment_rt,
                           young, old)
  ##CONSTRUCT METABRIC FEATURES DATA FRAME
  metabricData = data.frame(postMenopausal, stage, lymphNodesRemoved, NPI,
                            p53status_wt, p53status_mut, intClust_1,intClust_2,
                            intClust_3, intClust_4, intClust_5, intClust_6,
                            intClust_7, intClust_8, intClust_9, intClust_10,                        
                            group_1, group_2, group_3,group_4, group_other,                      
                            pam50_basal, pam50_her2, pam50_lumA, pam50_lumB, pam50_normal,
                            cellularity_low, cellularity_mod, cellularity_high,                        
                            site_1, site_2, site_3, site_4, site_5)
  ##OPTION TO ADD METABRIC-ONLY FEATURES
  if (is.null(site_1)==FALSE)
    outputData = cbind(outputData, metabricData)
  #SUBSET OF FEATURES
  outputData = data.frame(age, size, lymphNodes, treatment_ct, treatment_ht, treatment_rt,
                          pam50_lumA, grade, p53status_wt, her2.snp6.gain, postMenopausal,
                          prStatus.expr, young, old, pam50_lumB)
  ##PATCH MISSING DATA
  #   for (i in 1:ncol(outputData)){
  #     currentData        = outputData[,i]
  #     patch              = which(is.na(currentData))
  #     currentData[patch] = median(currentData, na.rm=TRUE)
  #     outputData[,i]     = currentData
  #   }
  return(outputData)
}

# train a RSF survival model
train.gene.rsf <- function(data, clinicalSurvData, numSplit=0, numTrees=200){
  data$time   <- clinicalSurvData[,1]
  data$status <- clinicalSurvData[,2]
  model <- rsf(Surv(time,status) ~ ., data=as.data.frame(data), ntree=numTrees, nsplit=numSplit)
  return(model)
}


predict.gene.rsf <- function(model, data){
  predictObject <- predict(model, as.data.frame(data))
  predictions   <- as.numeric(predictObject$mortality)
  return(predictions)
}

GetPredictCI <- function(pred, clinicalSurvData) {
  return(SurvivalModelPerformance$new(pred, clinicalSurvData)$getConcordanceIndex()$c.index)
}


##----------------------------------------------------------------------
## MODEL:  GLMNET ------------------------------------------------------
##----------------------------------------------------------------------
##METHOD TO TRAIN MODEL
## wrapper for glmnet with cox model from Rich Savage 
Train.glmnet <- function(data, targetValues){
  data                   <- as.data.frame(data)
  colnames(targetValues) <- c("time", "status")
  ##ADD AN OFFSET TERM
  nDataItems   <- nrow(data) 
  offsetVector <- 1 + numeric(nDataItems)
  data$offset  <- offsetVector
  ##TRAIN MODEL 1
  model.1           <- glmnet(as.matrix(data),    targetValues, family="cox", maxit=1e4, alpha=1)
  cv.fit.1          <- cv.glmnet(as.matrix(data), targetValues, family="cox", maxit=1e4, alpha=1)
  model.1$minLambda <- cv.fit.1$lambda.min
  ##TRAIN MODEL 2
  model.2           <- glmnet(as.matrix(data),    targetValues, family="cox", maxit=1e4, alpha=0.5)
  cv.fit.2          <- cv.glmnet(as.matrix(data), targetValues, family="cox", maxit=1e4, alpha=0.5)
  model.2$minLambda <- cv.fit.2$lambda.min
  ##TRAIN MODEL 3
  model.3 = Train.cox(data, targetValues)
  ##STORE MODELS TOGETHER AND RETURN
  model = list(model.1, model.2, model.3)
  return(model)
}

##METHOD TO MAKE PREDICTIONS USING THE MODEL
Predict.glmnet <- function(model, data){
  if (class(model)[1]!="uninitializedField"){
    data = as.data.frame(data)
    ##ADD AN OFFSET TERM
    nDataItems   <- nrow(data) 
    offsetVector <- 1 + numeric(nDataItems)
    data$offset  <- offsetVector
    ##MAKE PREDICTIONS                         
    ##print("Converting data.cl to a numerical matrix here may be an issue....")
    predictions.1 <- predict(model[[1]], as.matrix(data), type="response", s=model[[1]]$minLambda)
    predictions.2 <- predict(model[[2]], as.matrix(data), type="response", s=model[[2]]$minLambda)
    predictions.3 <- Predict.cox(model[[3]], data)
    predictions   <- as.numeric(predictions.1) + as.numeric(predictions.2)
    predictions   <- scale(predictions) + predictions.3
    predictions   <- scale(predictions)
  }
  else predictions <- numeric(nrow(data))
  return(predictions)
}

##----------------------------------------------------------------------
## MODEL:  K-NEAREST-NEIGHBOURS REGRESSION -----------------------------
##----------------------------------------------------------------------
## wrapper for kNN from Rich Savage
Predict.knn <- function(data.train, data.test, trainingTargets){
  ##FIND USEFUL VALUES  
  nDataItems  <- nrow(data.train)
  index       <- 1:nDataItems
  working     <- as.numeric(trainingTargets)
  timeToEvent <- working[index]
  status      <- working[-index]
  ##MAKE PREDICTIONS
  model       = knn.reg(as.matrix(data.train), as.matrix(data.test), timeToEvent, k=11)
  predictions = as.numeric(model$pred)
  predictions = -predictions
  predictions = scale(predictions)
  return(predictions)
}

PrintCI <- function(pred, clinicalSurvData, verbose=FALSE) {
  out <- NULL
  if (is.null(nrow(pred))) {
    #tmp <- SurvivalModelPerformance$new(pred, clinicalSurvData)$getConcordanceIndex()$c.index
    tmp <- concordance.index(x=pred, surv.time=clinicalSurvData[,1], surv.event=clinicalSurvData[,2], na.rm=TRUE, alpha= .05)$c.index
    
    if (verbose) cat(tmp,"\t", 1-tmp, "\n")
    out <- c(tmp, 1-tmp, max(tmp,1-tmp))
  } else if (nrow(pred)==nrow(clinicalSurvData)) {
    for (i in 1:ncol(pred)) {
      #tmp <- SurvivalModelPerformance$new(pred[,i], clinicalSurvData)$getConcordanceIndex()$c.index
      tmp <- concordance.index(x=pred[,i], surv.time=clinicalSurvData[,1], surv.event=clinicalSurvData[,2], na.rm=TRUE, alpha= .05)$c.index
      if (verbose) cat(i,"\t", tmp,"\t", 1-tmp, "\n")
      out <- rbind(out, c(i, tmp, 1-tmp,max(tmp,1-tmp)))
    }
  } else if (ncol(pred)==nrow(clinicalSurvData)) {
    for (i in 1:nrow(pred)) {
      #tmp <- SurvivalModelPerformance$new(pred[i,], clinicalSurvData)$getConcordanceIndex()$c.index
      tmp <- concordance.index(x=pred[i,], surv.time=clinicalSurvData[,1], surv.event=clinicalSurvData[,2], na.rm=TRUE, alpha= .05)$c.index
      if (verbose) cat(i,"\t", tmp,"\t", 1-tmp, "\n")    
      out <- rbind(out, c(i, tmp, 1-tmp,max(tmp,1-tmp)))
    }
  } else {
    cat("Error: sample size not equal\n")
  }
  return(out)
}


PrintExactCI <- function(pred, clinicalSurvData, verbose=FALSE) {
  out <- NULL
  if (is.null(nrow(pred))) {
    #tmp <- SurvivalModelPerformance$new(pred, clinicalSurvData)$getConcordanceIndex()$c.index
    tmp <- exactConcordanceIndex(pred, clinicalSurvData)$c.index
    
    if (verbose) cat(tmp,"\t", 1-tmp, "\n")
    out <- c(tmp, 1-tmp, max(tmp,1-tmp))
  } else if (nrow(pred)==nrow(clinicalSurvData)) {
    for (i in 1:ncol(pred)) {
      #tmp <- SurvivalModelPerformance$new(pred[,i], clinicalSurvData)$getConcordanceIndex()$c.index
      tmp <- exactConcordanceIndex(pred[,i], clinicalSurvData)$c.index
      if (verbose) cat(i,"\t", tmp,"\t", 1-tmp, "\n")
      out <- rbind(out, c(i, tmp, 1-tmp,max(tmp,1-tmp)))
    }
  } else if (ncol(pred)==nrow(clinicalSurvData)) {
    for (i in 1:nrow(pred)) {
      #tmp <- SurvivalModelPerformance$new(pred[i,], clinicalSurvData)$getConcordanceIndex()$c.index
      tmp <- exactConcordanceIndex(pred[i,], clinicalSurvData)$c.index
      if (verbose) cat(i,"\t", tmp,"\t", 1-tmp, "\n")    
      out <- rbind(out, c(i, tmp, 1-tmp,max(tmp,1-tmp)))
    }
  } else {
    cat("Error: sample size not equal\n")
  }
  return(out)
}


# CIBoost Class
CIBoost <- setRefClass(Class = "CIBoost",
                       fields = c("weights","gridWeight","maxIter","alpha","minIncr"),
                       methods = list(
                         initialize = function(maxIter=10,alpha=1,minIncr=0,gridMax=0.5,gridNum=50,...){
                           .self$gridWeight <- (ppoints(gridNum)-0.5)*2*gridMax  # from -gridMax to gridMax
                           .self$maxIter <- maxIter
                           .self$alpha <- alpha  # shrinkage parameter
                           .self$minIncr <- minIncr
                           .self$weights <- NULL
                           return(.self)                           
                         },
                         
                         # set up initial weight
                         initWeightsWithMaxCI = function(pred,survData) {
                           # find the featue with max CI
                           tmp <- NULL
                           for (i in 1:ncol(pred)) {
                             s <- SurvivalModelPerformance$new(pred[,i], survData)$getConcordanceIndex()$c.index
                             if (s>0.5) tmp <- rbind(tmp, c(i,1,s))
                             else tmp <- rbind(tmp, c(i,-1,1-s))
                           }
                           
                           currID <- which.max(tmp[,3])
                           currPred <- pred[,currID]
                           if (tmp[currID,2]==-1) currPred <- -currPred  # reverse the sign
                           
                           .self$weights <- tmp[currID,]  # init weights with the best pred                            
                         },
                         
                         # set up initial weight
                         initWeightsWithCox = function(pred,survData, useAIC=FALSE) {
                           # find the featue with max CI
                           
                           if (useAIC) {
                             upper <- terms(survData~(.), data=as.data.frame(pred))
                             fit <- step(coxph(survData~1, data=as.data.frame(pred)),scope=upper, direction="both", k=2, trace=FALSE)
                           } else {
                             fit <- coxph(survData~.,as.data.frame(pred))
                           }
                           coef <- fit$coef
                           
                           usedFeatures <- attr(fit$terms,"term.labels")
                           allFeatures <- colnames(pred)
                           
                           .self$weights <- NULL
                           for (f in usedFeatures) {
                             i <- which(allFeatures==f)
                             .self$weights <- rbind(.self$weights, c(i, fit$coef[f], 0))                             
                           }                  
                           
                           p <- customPredict(pred)
                           sig <- sd(p)
                           for (i in 1:nrow(.self$weights)) {
                             .self$weights[i,2] <- .self$weights[i,2]/sig
                           }
                         },
                         
                         # find weights based on CI 
                         customTrain = function(pred, survData, initUseAIC=FALSE, initWithMaxCI=FALSE) {
                           # pred: numSamples x numFeatures
                           
                           pred <- scale(pred)  # scale the column of each eacture
                           
                           # init weights
                           if (initWithMaxCI) {
                             initWeightsWithMaxCI(pred,survData)
                           } else {
                             initWeightsWithCox(pred,survData,initUseAIC)
                           }
                           
                           # calc currPred from init weights
                           currPred <- customPredict(pred)
                           bestCI <- SurvivalModelPerformance$new(currPred, survData)$getConcordanceIndex()$c.index
                           
                           n <- .self$maxIter   # remaining number of iterations
                           while (n>0) {
                             n <- n-1
                             
                             tmpList <- NULL
                             for (i in 1:ncol(pred)) {
                               w <- findWeight(currPred, pred[,i], survData)
                               newPred <- currPred + w*pred[,i]
                               tmpList <- rbind(tmpList, c(i,w,SurvivalModelPerformance$new(newPred, survData)$getConcordanceIndex()$c.index))
                             }
                             if (is.null(tmpList)) break
                             
                             id <- which.max(tmpList[,3])
                             # if increment is too small, break 
                             if ( tmpList[id,3] - bestCI < .self$minIncr) break
                             
                             bestCI <- tmpList[id,3]
                             .self$weights <- rbind(.self$weights, tmpList[id,])
                             currPred <- currPred + alpha*tmpList[id,2]*pred[,tmpList[id,1]]                                              
                           }                                   
                         },
                         
                         customPredict = function(pred) {
                           pred <- scale(pred)
                           p <- 0
                           if (is.null(nrow(.self$weights))) return(pred[,.self$weights[1]]*.self$weights[2])
                           
                           for (i in 1:nrow(.self$weights)) {
                             p <- p + pred[,.self$weights[i,1]] * .self$weights[i,2] * alpha
                             #p <- p/sd(p)
                           }
                           return(p)
                         },
                         
                         # find pairwise weight 
                         findWeight = function(predA, predB, survData) {
                           # assume both features are normalized to be mean 0 and var 1                           
                           tmp <- numeric(length(.self$gridWeight))
                           for (i in 1:length(tmp))    
                             tmp[i] <- SurvivalModelPerformance$new(predA + .self$gridWeight[i]*predB, survData)$getConcordanceIndex()$c.index 
                           
                           return(.self$gridWeight[which.max(tmp)])                                                        
                         }
                         
                       ) # end of list of methods
                       
) # end of class definition



# definition of ConcordanceIndex class
ConcordanceIndex <- setRefClass(Class = "ConcordanceIndex",
                                fields = c("time", "event", "ordered.index","rank", "event.index", "total.pairs", "size", "alpha"),
                                # ordered.index:  index of the surv orderd by surv time 
                                # rank: rank of surv samples based on surv time from low to high
                                # event.index: index of ordered.index that are events
                                methods = list(
                                  initialize = function(surv,...){
                                    .self$ordered.index <- order(surv[,1]) # ordered by time
                                    .self$size <- length(.self$ordered.index)
                                    
                                    .self$event <- as.numeric(surv[,2])
                                    .self$time <- as.numeric(surv[,1])
                                    
                                    # find the rank of each sample based on its surv time 
                                    .self$rank <- order(.self$ordered.index)
                                    
                                    .self$alpha <- 1  # slope of the sigmoidal function
                                    
                                    # find the index of ordered.index vector that that underlying sample has an event
                                    .self$event.index <- NULL
                                    .self$total.pairs <- 0
                                    for (i in 1:.self$size) {
                                      if (surv[.self$ordered.index[i],2]==1) {
                                        .self$event.index <- c(.self$event.index,i)
                                        .self$total.pairs <- .self$total.pairs + .self$size - i
                                      }                              
                                    }
                                    
                                    return(.self)                           
                                  },
                                  
                                  
                                  
                                  # f: a vector of risk predictions, higher risk -> shorter survival time
                                  # calculate the gradient of smoothed CI function
                                  gradient = function(f) {
                                    if (length(f) != .self$size) 
                                      stop("pred does not have the same size")
                                    
                                    # calculate exp function
                                    y <- exp(-.self$alpha * f)
                                    
                                    # calculate gradient
                                    g <- numeric(.self$size)
                                    for (k in 1:.self$size) {
                                      r <- .self$rank[k]  # corresponding location in ordered.index
                                      
                                      g[k] <- 0
                                      # if k is an event sample, find all (k,j) pairs
                                      if (.self$event[k]==1 & r<.self$size) {
                                        # find all samples at the right side of k: (k,j) pairs
                                        j <- .self$ordered.index[(r+1):.self$size]
                                        tmp <- y[k]/y[j]
                                        g[k] <- g[k] + sum(.self$alpha*tmp/(1+tmp)^2)
                                      }
                                      
                                      # find all (i,k) pairs with i being events
                                      if (r>1) {
                                        foo <- logical(.self$size)
                                        foo[.self$ordered.index[1:r-1]] <- TRUE
                                        
                                        i <- which(.self$event==1 & foo)  # both time < Ti and having an event
                                        if (length(i)>0) {
                                          tmp <- y[i]/y[k]
                                          g[k] <- g[k] + sum(-.self$alpha*tmp/(1+tmp)^2)
                                        }                                 
                                      }
                                      
                                    } # end of k loop
                                    g <- g/.self$total.pairs
                                    return(g)
                                  },
                                  
                                  # calcualte optimal step size based on risk prediction f and gradeint g
                                  calc.stepsize = function(f, g, max.stepsize=5, ci.resol=0.0001) {
                                    #delta <- c(0, ppoints(1/resol))*max.step
                                    
                                    left <- 0
                                    right <- max.stepsize
                                    ci.left <- get.smooth.ci(f+left*g)
                                    ci.right <- get.smooth.ci(f+right*g)
                                    
                                    while (abs(ci.left-ci.right) > ci.resol) {
                                      mid <- (left+right)/2
                                      ci <- get.smooth.ci(f+mid*g)
                                      if (ci.left < ci.right) {
                                        ci.left <- ci
                                        left <- mid
                                      } else {
                                        ci.right <- ci
                                        right <- mid
                                      }
                                      #cat(mid,"\t",ci,"\n")
                                    }
                                    
                                    if (ci.left > ci.right) 
                                      return(left)
                                    else
                                      return(right)                      
                                  },
                                  
                                  
                                  
                                  # f: a vector of risk predictions, higher risk -> shorter survival time
                                  # calculate the smoothed CI function
                                  get.smooth.ci = function(f) {
                                    if (length(f) != .self$size) 
                                      stop("pred does not have the same size")
                                    
                                    # calculate exp function
                                    y <- exp(-.self$alpha * f)
                                    
                                    sci <- 0
                                    for (k in .self$event.index) {
                                      if (k==.self$size) break
                                      i <- .self$ordered.index[k]               
                                      j <- .self$ordered.index[(k+1):.self$size]
                                      tmp <- y[i]/y[j]
                                      sci <- sci + sum(1/(1+tmp))    
                                    }
                                    
                                    return(sci/.self$total.pairs)
                                  },
                                  
                                  
                                  # calculate CI given a prediction
                                  # pred: a vector of risk predictions, higher risk -> shorter survival time
                                  get.ci = function(pred) {
                                    if (length(pred) != .self$size) 
                                      stop("pred does not have the same size")
                                    
                                    ci <- concordance.index(x=pred, surv.time=.self$time, surv.event=.self$event, na.rm=TRUE, alpha= .05)$c.index
                                    return(ci)             
                                  },
                                  
                                  # calculate CI given a prediction
                                  # pred: a vector of risk predictions, higher risk -> shorter survival time
                                  get.ci.my = function(pred) {
                                    if (length(pred) != .self$size) 
                                      stop("pred does not have the same size")
                                    
                                    n <- 0
                                    for (k in .self$event.index) {
                                      i <- .self$ordered.index[k]
                                      if (k<.self$size) {
                                        j <- .self$ordered.index[(k+1):.self$size]
                                        n <- n + sum(pred[j]<pred[i])  
                                      }
                                    }
                                    
                                    return(n/.self$total.pairs)             
                                  }
                                  
                                ) # end of list of methods
                                
) # end of class definition


# definition of ConcordanceIndexBoost class
ConcordanceIndexBoost <- setRefClass(Class = "ConcordanceIndexBoost",
                                     fields = c("ntrees","shrinkage","maxdepth","trees", "model.cox","max.stepsize","alpha"),
                                     methods = list(
                                       initialize = function(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1,...){                                                     
                                         .self$ntrees <- ntrees
                                         .self$shrinkage <- shrinkage  # shrinkage paramters
                                         .self$maxdepth <- maxdepth
                                         .self$max.stepsize <- max.stepsize
                                         .self$alpha <- alpha
                                         .self$trees <- NULL
                                         
                                         return(.self)                           
                                       },
                                       
                                       customTrain = function(survData, X) {
                                         # initialize with cox prediction
                                         upper <- terms(survData~(.), data = X)
                                         .self$model.cox <- step(coxph(survData~1, data=X),scope=upper, direction="both", k=2, trace=FALSE)
                                         p <- predict(.self$model.cox)
                                         p <- (p-mean(p))/sd(p)
                                         
                                         # init CI object with surv data  
                                         ci <- ConcordanceIndex$new(survData) 
                                         ci$alpha <- .self$alpha
                                         
                                         curr.sci <- ci$get.smooth.ci(p)
                                         
                                         .self$trees <- NULL
                                         for (iter in 1:ntrees) {
                                           g <- ci$gradient(p)
                                           delta <- ci$calc.stepsize(p, g, max.stepsize=.self$max.stepsize, ci.resol=0.0001)
                                           
                                           Y <- .self$shrinkage * delta * g
                                           
                                           fit <- rpart(Y~., data=X, method="anova",control=rpart.control(cp=0.001,minsplit=20,maxdepth=.self$maxdepth))
                                           
                                           .self$trees <- c(.self$trees, list(fit))
                                           p <- p + predict(fit)
                                           
                                           sci <- ci$get.smooth.ci(p)
                                           if (abs(sci-curr.sci) < 0.00001) break
                                           
                                           cat(iter,"\t", delta,"\t", sci,"\t", ci$get.ci(p),"\n")
                                           curr.sci <- sci
                                         }
                                       },
                                       
                                       
                                       
                                       customPredict = function(X) {
                                         p <- predict(.self$model.cox,X)
                                         p <- (p-mean(p))/sd(p)
                                         
                                         for (i in 1:length(trees)) {
                                           p <- p + predict(.self$trees[[i]], X)                                  
                                         }
                                         return(p)
                                       }
                                       
                                     ) # end of list of methods
                                     
) # end of class definition


###################################################################################
# using only clinical features 
EnsembleModelClncFeatures <- setRefClass(Class  = "EnsembleModelClncFeatures",
                                         contains = "PredictiveModel",
                                         fields   = c("model", "predictions", "use_osloval", "submit.model.index"),
                                         methods  = list(
                                           
                                           initialize = function(...) {
                                             .self$use_osloval <- TRUE
                                             .self$submit.model.index <- 4
                                             return(.self)
                                           },
                                           
                                           customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                             if(class(clinicalSurvData) != "Surv") {
                                               stop("Expecting 'responseData' object of type 'Surv'")
                                             } 
                                             rm(copyData)
                                             rm(exprData)
                                             
                                             clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                          
                                             clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                             
                                             
                                             X <- clinical
                                             upper <- terms(clinicalSurvData~(.), data = X)
                                             
                                             # AIC cox model
                                             cat("Train AIC Cox model ..."); flush.console()
                                             mcox <- step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=2, trace=FALSE)
                                             cat("done!\n");flush.console()
                                             
                                             X <- X[, attr(mcox$term, "term.labels")]  # use only cox AIC selected features
                                             
                                             
                                             # gbm model using cox AIC selected features
                                             cat("Train gbm model..."); flush.console()                                             
                                             mgbm <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=10, cv.folds=5, verbose=F, seed=53) 
                                             cat("done!\n");flush.console()
                                             
                                             for (i in 1:ncol(X)) {
                                               X[,i] <- as.numeric(X[,i])  # convert X to numeric values
                                             }
                                             
                                             # ci model
                                             cat("Train ci model ..."); flush.console()
                                             mci <- CIBoost$new(maxIter=10,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                             #mci$customTrain(X, clinicalSurvData, initUseAIC=TRUE)
                                             mci$customTrain(X, clinicalSurvData)
                                             cat("Done\n")
                                             
                                             # cib model
                                             cat("Train cib model ..."); flush.console()
                                             mcib <- ConcordanceIndexBoost$new(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1)
                                             mcib$customTrain(clinicalSurvData, X)
                                             cat("Done\n")                 
                                             
                                             .self$model <- list(mcox=mcox, mgbm=mgbm, mcib=mcib, mci=mci)
                                           },
                                           
                                           
                                           
                                           customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                             rm(copyData)
                                             rm(exprData)
                                             
                                             clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                            
                                             clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)                                                
                                             
                                             # cox AIC model prediction 
                                             X <- clinical
                                             p1 <- predict(.self$model$mcox, X)
                                             
                                             X <- X[, attr(.self$model$mcox$term, "term.labels")]
                                             
                                             # gbm model prediction                                                                                     
                                             best.iter <- gbm.perf(.self$model$mgbm, method="cv", plot.it=FALSE)
                                             p2 <- predict.gbm(.self$model$mgbm, X, best.iter)
                                             
                                             for (i in 1:ncol(X)) {
                                               X[,i] <- as.numeric(X[,i])  # convert X to numeric values
                                             }
                                             
                                             # ci model prediction
                                             p3 <- .self$model$mci$customPredict(X)                      
                                             
                                             # cib model prediction
                                             p4 <- .self$model$mcib$customPredict(X)                      
                                             
                                             pz <- cbind(p1,p2,p3,p4)
                                             
                                             .self$predictions <- pz
                                             
                                             pz <- scale(pz)
                                             #p <- apply(pz, 1, mean) # combine two models
                                             
                                             
                                             if (length(.self$submit.model.index)==1) 
                                               p <- pz[,.self$submit.model.index]
                                             else 
                                               p <- apply(pz[,.self$submit.model.index],1,mean)
                                             
                                             
                                             return (p)
                                           }
                                           
                                         )
)



###########################

EnsembleModelClncGeneFeatures <- setRefClass(Class  = "EnsembleModelClncGeneFeatures",
                                             contains = "PredictiveModel",
                                             fields   = c("model", "predictions", "selectedGenes", "use_osloval", "submit.model.index","combined.features"),
                                             methods  = list(
                                               
                                               initialize = function(...) {                                            
                                                 .self$use_osloval <- TRUE
                                                 .self$submit.model.index <- c(3,4)                                                  
                                                 .self$selectedGenes <- c('ILMN_1683450', 'ILMN_2392472', 'ILMN_1700337', 'ILMN_1781943', 'ILMN_1772686', 
                                                                          'ILMN_2077550', 'ILMN_1796949', 'ILMN_1716279', 'ILMN_1808071', 'ILMN_1714730',
                                                                          'ILMN_1785570', 'ILMN_1731070', 'ILMN_1801257', 'ILMN_1736176', 'ILMN_1763907')
                                                 return(.self)
                                               },
                                               
                                               
                                               customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                                 if(class(clinicalSurvData) != "Surv") {
                                                   stop("Expecting 'responseData' object of type 'Surv'")
                                                 } 
                                                 rm(copyData)
                                                 
                                                 clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                          
                                                 clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                                 
                                                 # extract expression of selectedGenes
                                                 selectedExprData <- ExtractGeneExprData(exprData, .self$selectedGenes) 
                                                 
                                                 
                                                 # clnc AIC cox model
                                                 X <- clinical
                                                 upper <- terms(clinicalSurvData~(.), data = X)
                                                 clnc.cox <- step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=3, trace=FALSE)
                                                 
                                                 # gene AIC cox model
                                                 X <- cbind(clinical[,attr(clnc.cox$terms, "term.labels")], selectedExprData)
                                                 upper <- terms(clinicalSurvData~(.), data = X)
                                                 expr.cox <- step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=1, trace=FALSE)
                                                 
                                                 #.self$combined.features <- c(attr(clnc.cox$terms, "term.labels"), attr(expr.cox$terms,"term.labels"))
                                                 .self$combined.features <- attr(expr.cox$terms,"term.labels")
                                                 
                                                 # joined features
                                                 X = cbind(selectedExprData, clinical)
                                                 X = X[, .self$combined.features]
                                                 
                                                 # train cox AIC model
                                                 upper <- terms(clinicalSurvData~(.), data = X)
                                                 mcox <- step(coxph(clinicalSurvData~., data=X), scope=upper, direction="both", k=2, trace=FALSE)
                                                 
                                                 # train gbm model
                                                 mgbm <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=6, cv.folds=5, verbose=F, seed=123) 
                                                 
                                                 
                                                 for (i in 1:ncol(X)) {
                                                   X[,i] <- as.numeric(X[,i])
                                                 }    
                                                 
                                                 # ci model
                                                 cat("Train ci model ..."); flush.console()
                                                 mci <- CIBoost$new(maxIter=5,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                                 #mci$customTrain(X, clinicalSurvData,initWithMaxCI=TRUE)  #init with max CI genes
                                                 mci$customTrain(X, clinicalSurvData,initUseAIC=TRUE)  
                                                 cat("Done\n")
                                                 
                                                 # cib model
                                                 cat("Train cib model ..."); flush.console()
                                                 mcib <- ConcordanceIndexBoost$new(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1)
                                                 mcib$customTrain(clinicalSurvData, X)
                                                 cat("Done\n")                 
                                                 
                                                 .self$model <- list(mcox=mcox, mgbm=mgbm, mcib=mcib, mci=mci)
                                               },
                                               
                                               
                                               
                                               customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                                 rm(copyData)
                                                 
                                                 clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                          
                                                 clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                                 
                                                 # extract expression of selectedGenes
                                                 selectedExprData <- ExtractGeneExprData(exprData, .self$selectedGenes) 
                                                 
                                                 # joined features
                                                 X = cbind(selectedExprData, clinical)                                
                                                 X = X[, .self$combined.features]
                                                 
                                                 #X = X[, attr(.self$model$mcox$term, "term.labels")]
                                                 
                                                 # cox AIC model prediction                 
                                                 p1 <- predict(.self$model$mcox, X)
                                                 
                                                 # gbm model prediction                                                                                     
                                                 best.iter <- gbm.perf(.self$model$mgbm, method="cv", plot.it=FALSE)
                                                 p2 <- predict.gbm(.self$model$mgbm, X, best.iter)
                                                 
                                                 
                                                 for (i in 1:ncol(X)) {
                                                   X[,i] <- as.numeric(X[,i])
                                                 }    
                                                 
                                                 # ci model prediction
                                                 p3 <- .self$model$mci$customPredict(X)                                                                  
                                                 
                                                 # cib model prediction
                                                 p4 <- .self$model$mcib$customPredict(X)                      
                                                 
                                                 pz <- cbind(p1,p2,p3,p4)                                                                                         
                                                 pz <- scale(pz)
                                                 .self$predictions <- pz
                                                 
                                                 if (length(.self$submit.model.index)==1) 
                                                   p <- pz[,.self$submit.model.index]
                                                 else 
                                                   p <- apply(pz[,.self$submit.model.index],1,mean)
                                                 
                                                 
                                                 return (p)
                                               }
                                               
                                             )
)

#########################################################################
EnsembleModelGeneFeatures <- setRefClass(Class  = "EnsembleModelGeneFeatures",
                                         contains = "PredictiveModel",
                                         fields   = c("model", "predictions", "selectedGenes", "use_osloval", "submit.model.index"),
                                         methods  = list(
                                           
                                           initialize = function(...) {                                            
                                             .self$use_osloval <- TRUE
                                             .self$submit.model.index <- 3                                                  
                                             .self$selectedGenes <- c('ILMN_1683450', 'ILMN_2392472', 'ILMN_1700337', 'ILMN_1781943', 'ILMN_1772686', 
                                                                      'ILMN_2077550', 'ILMN_1796949', 'ILMN_1716279', 'ILMN_1808071', 'ILMN_1714730',
                                                                      'ILMN_1785570', 'ILMN_1731070', 'ILMN_1801257', 'ILMN_1736176', 'ILMN_1763907')
                                             #.self$selectedGenes <- .self$selectedGenes[1:10]
                                             return(.self)
                                           },
                                           
                                           
                                           customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                             if(class(clinicalSurvData) != "Surv") {
                                               stop("Expecting 'responseData' object of type 'Surv'")
                                             } 
                                             rm(copyData)
                                             
                                             # extract expression of selectedGenes
                                             X <- ExtractGeneExprData(exprData, .self$selectedGenes) 
                                             
                                             # train cox AIC model
                                             upper <- terms(clinicalSurvData~(.), data = X)
                                             mcox <- step(coxph(clinicalSurvData~., data=X), scope=upper, direction="both", k=2, trace=FALSE)
                                             
                                             # train gbm model
                                             mgbm <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=6, cv.folds=5, verbose=F, seed=123) 
                                             
                                             
                                             # ci model
                                             cat("Train ci model ..."); flush.console()
                                             mci <- CIBoost$new(maxIter=4,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                             mci$customTrain(X, clinicalSurvData,initWithMaxCI=TRUE)  #init with max CI genes
                                             cat("Done\n")
                                             
                                             label.ci <- colnames(X)[mci$weights[,1]]
                                             X <- X[, union(attr(mcox$term, "term.labels"),label.ci)]
                                             # cib model
                                             cat("Train cib model ..."); flush.console()
                                             mcib <- ConcordanceIndexBoost$new(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1)
                                             mcib$customTrain(clinicalSurvData, X)
                                             cat("Done\n")                 
                                             
                                             .self$model <- list(mcox=mcox, mgbm=mgbm, mcib=mcib, mci=mci)
                                           },
                                           
                                           
                                           
                                           customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                             rm(copyData)
                                             
                                             # extract expression of selectedGenes
                                             X <- ExtractGeneExprData(exprData, .self$selectedGenes) 
                                             
                                             #X = X[, attr(.self$model$mcox$term, "term.labels")]
                                             
                                             # cox AIC model prediction                 
                                             p1 <- predict(.self$model$mcox, X)
                                             
                                             # gbm model prediction                                                                                     
                                             best.iter <- gbm.perf(.self$model$mgbm, method="cv", plot.it=FALSE)
                                             p2 <- predict.gbm(.self$model$mgbm, X, best.iter)
                                             
                                             
                                             # ci model prediction
                                             p3 <- .self$model$mci$customPredict(X)                      
                                             
                                             label.ci <- colnames(X)[.self$model$mci$weights[,1]]
                                             X <- X[, union(attr(.self$model$mcox$term, "term.labels"),label.ci)]
                                             
                                             #X <- X[, attr(.self$model$mcox$term, "term.labels")]
                                             
                                             # cib model prediction
                                             p4 <- .self$model$mcib$customPredict(X)                      
                                             
                                             pz <- cbind(p1,p2,p3,p4)                                                                                         
                                             pz <- scale(pz)
                                             .self$predictions <- pz
                                             
                                             if (length(.self$submit.model.index)==1) 
                                               p <- pz[,.self$submit.model.index]
                                             else 
                                               p <- apply(pz[,.self$submit.model.index],1,mean)
                                             
                                             
                                             return (p)
                                           }
                                           
                                         )
)

################################################################################

# using metagene features from the Goldi model from the DreamBox package
EnsembleModelMetageneFeatures <- setRefClass(Class  = "EnsembleModelMetageneFeatures",
                                             contains = "PredictiveModel",
                                             fields   = c("model", "predictions", "attractome", "annot", "chosenProbes", "use_osloval", "submit.model.index"),
                                             methods  = list(
                                               
                                               initialize = function(...) {
                                                 data(attractome.minimalist)
                                                 .self$attractome <- attractome.minimalist
                                                 data(map)
                                                 .self$annot <- map                                                 
                                                 .self$use_osloval <- FALSE
                                                 .self$submit.model.index <- c(1,2)
                                                 return(.self)
                                               },
                                               
                                               customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                                 if(class(clinicalSurvData) != "Surv") {
                                                   stop("Expecting 'responseData' object of type 'Surv'")
                                                 } 
                                                 rm(copyData)
                                                 
                                                 clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                          
                                                 clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                                 
                                                 # create metagene space    
                                                 o = CreateMetageneSpace(exprs(exprData), .self$attractome, .self$annot)
                                                 meta = o$metaSpace
                                                 .self$chosenProbes = o$pbs
                                                 ls = meta["ls",]
                                                 meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                                 
                                                 #===== Conditioning mesenchymal transition metagene by lymph node status
                                                 idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                                 mes.lymphneg = meta["mt",] * idx
                                                 mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                                 meta = rbind(meta, mes.lymphneg)
                                                 
                                                 idx = (clinical[,"lymph_nodes_positive"] > 3)
                                                 ls.lymphpos = ls * idx
                                                 ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                                 meta = rbind(meta, ls.lymphpos)
                                                 
                                                 idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                                 ls.erneg = ls * idx
                                                 ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                                 meta = rbind(meta, ls.erneg)
                                                 
                                                 
                                                 X <- data.frame(t(meta))
                                                 
                                                 # AIC cox model
                                                 cat("Train AIC Cox model ..."); flush.console()
                                                 upper <- terms(clinicalSurvData~(.), data = X)
                                                 mcox <- step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=2, trace=FALSE)
                                                 cat("done!\n");flush.console()
                                                 
                                                 
                                                 
                                                 # gbm model 
                                                 cat("Train gbm model..."); flush.console()                                             
                                                 mgbm <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=6, cv.folds=5, verbose=F, seed=123) 
                                                 cat("done!\n");flush.console()
                                                 
                                                 #X <- X[, attr(mcox$term, "term.labels")]  # use only cox AIC selected features
                                                 
                                                 # ci model
                                                 cat("Train ci model ..."); flush.console()
                                                 mci <- CIBoost$new(maxIter=10,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                                 #mci$customTrain(X, clinicalSurvData, initUseAIC=TRUE)
                                                 mci$customTrain(X, clinicalSurvData)
                                                 cat("Done\n")
                                                 
                                                 # cib model
                                                 cat("Train cib model ..."); flush.console()
                                                 mcib <- ConcordanceIndexBoost$new(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1)
                                                 mcib$customTrain(clinicalSurvData, X)
                                                 cat("Done\n")                 
                                                 
                                                 .self$model <- list(mcox=mcox, mgbm=mgbm, mcib=mcib, mci=mci)
                                               },
                                               
                                               
                                               
                                               customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                                 rm(copyData)
                                                 
                                                 clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                            
                                                 clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)                                                
                                                 
                                                 meta = CreateMetageneSpace(exprs(exprData), chosenProbes = .self$chosenProbes)
                                                 ls = meta["ls",]
                                                 meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                                 #===== Conditioning mesenchymal transition metagene by lymph node status
                                                 idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                                 mes.lymphneg = meta["mt",] * idx
                                                 mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                                 meta = rbind(meta, mes.lymphneg)
                                                 
                                                 idx = (clinical[,"lymph_nodes_positive"] > 3)
                                                 ls.lymphpos = ls * idx
                                                 ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                                 meta = rbind(meta, ls.lymphpos)
                                                 
                                                 idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                                 ls.erneg = ls * idx
                                                 ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                                 meta = rbind(meta, ls.erneg)
                                                 
                                                 X <- data.frame(t(meta))
                                                 
                                                 # cox AIC model prediction                 
                                                 p1 <- predict(.self$model$mcox, X)
                                                 
                                                 
                                                 # gbm model prediction                                                                                     
                                                 best.iter <- gbm.perf(.self$model$mgbm, method="cv", plot.it=FALSE)
                                                 p2 <- predict.gbm(.self$model$mgbm, X, best.iter)
                                                 
                                                 
                                                 #X <- X[, attr(.self$model$mcox$term, "term.labels")]
                                                 
                                                 # ci model prediction
                                                 p3 <- .self$model$mci$customPredict(X)                      
                                                 
                                                 # cib model prediction
                                                 p4 <- .self$model$mcib$customPredict(X)                      
                                                 
                                                 pz <- cbind(p1,p2,p3,p4)
                                                 
                                                 .self$predictions <- pz
                                                 
                                                 pz <- scale(pz)
                                                 
                                                 if (length(.self$submit.model.index)==1) 
                                                   p <- pz[,.self$submit.model.index]
                                                 else 
                                                   p <- apply(pz[,.self$submit.model.index],1,mean)
                                                 
                                                 
                                                 return (p)
                                               }
                                               
                                             )
)

##########################################################################################################

# predict using nearest neighbor combing clnc and meta gene features 
# adapted from the Goldi model authored by Wei-Yi Cheng
EnsembleModelMetageneNN <- setRefClass(Class  = "EnsembleModelMetageneNN",
                                       contains = "PredictiveModel",
                                       fields   = c("model", "predictions", "attractome", "annot", "chosenProbes", "use_osloval"),
                                       methods  = list(
                                         
                                         initialize = function(...) {
                                           data(attractome.minimalist)
                                           .self$attractome <- attractome.minimalist
                                           data(map)
                                           .self$annot <- map                                                 
                                           .self$use_osloval <- TRUE
                                           return(.self)
                                         },
                                         
                                         customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                           if(class(clinicalSurvData) != "Surv") {
                                             stop("Expecting 'responseData' object of type 'Surv'")
                                           } 
                                           rm(copyData)
                                           
                                           clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)  
                                           
                                           if (.self$use_osloval==FALSE) {
                                             used_col <- NULL
                                             for (c in colnames(clnc)) {
                                               if (substr(c,1,3)!='NOT') used_col <- c(used_col, c)
                                             }
                                             clnc <- clnc[,used_col]
                                           }
                                           
                                           clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                           
                                           # create metagene space    
                                           o = CreateMetageneSpace(exprs(exprData), .self$attractome, .self$annot)
                                           meta = o$metaSpace
                                           .self$chosenProbes = o$pbs
                                           ls = meta["ls",]
                                           meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                           
                                           #===== Conditioning mesenchymal transition metagene by lymph node status
                                           idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                           mes.lymphneg = meta["mt",] * idx
                                           mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                           meta = rbind(meta, mes.lymphneg)
                                           
                                           idx = (clinical[,"lymph_nodes_positive"] > 3)
                                           ls.lymphpos = ls * idx
                                           ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                           meta = rbind(meta, ls.lymphpos)
                                           
                                           idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                           ls.erneg = ls * idx
                                           ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                           meta = rbind(meta, ls.erneg)
                                           
                                           
                                           
                                           t = clinicalSurvData[,1]
                                           defSurvSamples = which(clinicalSurvData[,2]==1 | clinicalSurvData[,1] > 365 * 10)
                                           ccdi = getAllCCDIWz(meta, clinicalSurvData)
                                           idx = ccdi[c("er", "mitotic", "puf60", "erbb2", "chr7p11.2", "ls", "mt", "chr15q26.1")]
                                           
                                           knnmodel = list()
                                           knnmodel$x.train = list(meta=meta[names(idx), defSurvSamples], time=t[defSurvSamples], concordance=idx)
                                           
                                           knnmodel$c.train = preproClncKNN(clnc, clinicalSurvData, ccdi.upper=0.6, ccdi.lower=0.4)
                                           
                                           .self$model <- list(knnmodel=knnmodel)
                                         },
                                         
                                         
                                         
                                         customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                           rm(copyData)
                                           
                                           clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                            
                                           clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)                                                
                                           
                                           meta = CreateMetageneSpace(exprs(exprData), chosenProbes = .self$chosenProbes)
                                           ls = meta["ls",]
                                           meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                           #===== Conditioning mesenchymal transition metagene by lymph node status
                                           idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                           mes.lymphneg = meta["mt",] * idx
                                           mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                           meta = rbind(meta, mes.lymphneg)
                                           
                                           idx = (clinical[,"lymph_nodes_positive"] > 3)
                                           ls.lymphpos = ls * idx
                                           ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                           meta = rbind(meta, ls.lymphpos)
                                           
                                           idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                           ls.erneg = ls * idx
                                           ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                           meta = rbind(meta, ls.erneg)
                                           
                                           # Predicting using KNN model 
                                           knnmodel = .self$model$knnmodel
                                           qX = meta[names(knnmodel$x.train$concordance),]
                                           qC = clnc[,colnames(knnmodel$c.train$clinical)]
                                           qC = t(preproClncKNN(qC, isFactorIn=knnmodel$c.train$isFactor, dwIn=knnmodel$c.train$distWeight)$clinical)
                                           wvec = c(abs(knnmodel$x.train$concordance-0.5), abs(knnmodel$c.train$concordance-0.5))
                                           qAll = rbind(qX, qC)
                                           trainDB = rbind(knnmodel$x.train$meta, t(knnmodel$c.train$clinical))
                                           trainTime = knnmodel$x.train$time
                                           out=ewknn.predict(trainDB, trainTime, qAll, wvec, k=floor(0.1*ncol(trainDB)))
                                           p = 365/out
                                           .self$predictions <- p
                                           
                                           return (p)
                                         }
                                         
                                       )
)


#################################################################################################

# using the minimalist features from the Goldi model 
# train using CIBoost and ConcordanceIndexBoost
EnsembleModelMinFeatures <- setRefClass(Class  = "EnsembleModelMinFeatures",
                                        contains = "PredictiveModel",
                                        fields   = c("model", "attractome", "annot", "predictions", "chosenProbes", "use_osloval"),
                                        methods  = list(
                                          
                                          initialize = function(...) {
                                            data(attractome.minimalist)
                                            .self$attractome <- attractome.minimalist
                                            data(map)
                                            .self$annot <- map
                                            .self$use_osloval <- TRUE
                                            return(.self)
                                          },
                                          
                                          customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                            if(class(clinicalSurvData) != "Surv") {
                                              stop("Expecting 'responseData' object of type 'Surv'")
                                            } 
                                            rm(copyData)
                                            
                                            clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                          
                                            clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                            
                                            # create metagene space    
                                            o = CreateMetageneSpace(exprs(exprData), .self$attractome, .self$annot)
                                            meta = o$metaSpace
                                            .self$chosenProbes = o$pbs
                                            ls = meta["ls",]
                                            meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                            
                                            #===== Conditioning mesenchymal transition metagene by lymph node status
                                            idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                            mes.lymphneg = meta["mt",] * idx
                                            mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                            meta = rbind(meta, mes.lymphneg)
                                            
                                            idx = (clinical[,"lymph_nodes_positive"] > 3)
                                            ls.lymphpos = ls * idx
                                            ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                            meta = rbind(meta, ls.lymphpos)
                                            
                                            idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                            ls.erneg = ls * idx
                                            ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                            meta = rbind(meta, ls.erneg)
                                            
                                            rm(exprData)
                                            
                                            
                                            # minimalist features from the Goldi model
                                            X = data.frame( cbind(meta["mitotic",], meta["ls.erneg",], clinical$lymph_nodes_positive, meta["mes.lymphneg",], meta["susd3",], clinical$age_at_diagnosis, clinical$tr.RT, clinical$tr.CT) )
                                            colnames(X) = c("CIN", "LYM_ERNeg", "lymNum", "MES_lymNumNeg", "SUSD3", "age", "tr.RT", "tr.CT")
                                            
                                            # ci model
                                            cat("Train model.ci ..."); flush.console()
                                            model.ci <- CIBoost$new(maxIter=10,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                            model.ci$customTrain(X, clinicalSurvData)
                                            cat("Done\n")
                                            
                                            # cib model
                                            cat("Train model.cib ..."); flush.console()
                                            model.cib <- ConcordanceIndexBoost$new(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1)
                                            model.cib$customTrain(clinicalSurvData, X)
                                            cat("Done\n")                 
                                            
                                            .self$model <- list(model.cib=model.cib, model.ci=model.ci)
                                          },
                                          
                                          
                                          
                                          customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                            
                                            rm(copyData)
                                            
                                            clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                            
                                            clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                            
                                            meta = CreateMetageneSpace(exprs(exprData), chosenProbes = .self$chosenProbes)
                                            ls = meta["ls",]
                                            meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                            #===== Conditioning mesenchymal transition metagene by lymph node status
                                            idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                            mes.lymphneg = meta["mt",] * idx
                                            mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                            meta = rbind(meta, mes.lymphneg)
                                            
                                            idx = (clinical[,"lymph_nodes_positive"] > 3)
                                            ls.lymphpos = ls * idx
                                            ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                            meta = rbind(meta, ls.lymphpos)
                                            
                                            idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                            ls.erneg = ls * idx
                                            ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                            meta = rbind(meta, ls.erneg)
                                            
                                            rm(exprData)
                                            
                                            # predictions
                                            X = data.frame( cbind(meta["mitotic",], meta["ls.erneg",], clinical$lymph_nodes_positive, meta["mes.lymphneg",], meta["susd3",], clinical$age_at_diagnosis, clinical$tr.RT, clinical$tr.CT, clinical$size) )
                                            colnames(X) = c("CIN", "LYM_ERNeg", "lymNum", "MES_lymNumNeg", "SUSD3", "age", "tr.RT", "tr.CT", "size")
                                            
                                            
                                            p1 <- .self$model$model.ci$customPredict(X)                      
                                            p2 <- .self$model$model.cib$customPredict(X)                      
                                            pz <- cbind(p1,p2)
                                            
                                            .self$predictions <- pz
                                            
                                            pz <- scale(pz)
                                            p <- apply(pz, 1, mean) # combine two models
                                            return (p)
                                          }
                                          
                                        )
)


###########################################################################################################################

# using left over features from cox AIC clnc and expr models, adapted from the Goldi model
# train using CIBoost and ConcordanceIndexBoost
EnsembleModelSecondaryFeatures <- setRefClass(Class  = "EnsembleModelSecondaryFeatures",
                                              contains = "PredictiveModel",
                                              fields   = c("model", "predictions", "attractome", "annot", "chosenProbes", "use_osloval", "submit.model.index"),
                                              methods  = list(
                                                
                                                initialize = function(...) {
                                                  data(attractome.minimalist)
                                                  .self$attractome <- attractome.minimalist
                                                  data(map)
                                                  .self$annot <- map                                                 
                                                  .self$use_osloval <- TRUE
                                                  .self$submit.model.index <- c(2,3)
                                                  return(.self)
                                                },
                                                
                                                customTrain = function(exprData, copyData, clinicalFeaturesData,clinicalSurvData,...) {
                                                  if(class(clinicalSurvData) != "Surv") {
                                                    stop("Expecting 'responseData' object of type 'Surv'")
                                                  } 
                                                  rm(copyData)
                                                  
                                                  clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                          
                                                  clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)   
                                                  
                                                  # create metagene space    
                                                  o = CreateMetageneSpace(exprs(exprData), .self$attractome, .self$annot)
                                                  meta = o$metaSpace
                                                  .self$chosenProbes = o$pbs
                                                  ls = meta["ls",]
                                                  meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                                  
                                                  #===== Conditioning mesenchymal transition metagene by lymph node status
                                                  idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                                  mes.lymphneg = meta["mt",] * idx
                                                  mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                                  meta = rbind(meta, mes.lymphneg)
                                                  
                                                  idx = (clinical[,"lymph_nodes_positive"] > 3)
                                                  ls.lymphpos = ls * idx
                                                  ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                                  meta = rbind(meta, ls.lymphpos)
                                                  
                                                  idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                                  ls.erneg = ls * idx
                                                  ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                                  meta = rbind(meta, ls.erneg)
                                                  
                                                  
                                                  
                                                  # clnc AIC cox model
                                                  X <- clinical
                                                  upper <- terms(clinicalSurvData~(.), data = X)
                                                  clnc.cox <- step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=2, trace=FALSE)
                                                  
                                                  # gene AIC cox model
                                                  X <- data.frame(t(meta))
                                                  upper <- terms(clinicalSurvData~(.), data = X)
                                                  expr.cox <- step(coxph(clinicalSurvData~1, data=X),scope=upper, direction="both", k=2, trace=FALSE)
                                                  
                                                  # remaining features
                                                  X = cbind(data.frame(t(meta)), clinical)
                                                  X = X[, setdiff(colnames(X), c(  attr(clnc.cox$terms, "term.labels"), attr(expr.cox$terms,"term.labels")   ))]
                                                  upper = terms(clinicalSurvData~(.), data = X)
                                                  mcox = step(coxph(clinicalSurvData~., data=X), scope=upper, direction="both", k=2, trace=FALSE)
                                                  
                                                  X = X[, attr(mcox$term, "term.labels")]
                                                  mgbm <- gbm.cvrun(clinicalSurvData~., data=X ,distribution="coxph", shrinkage=0.002, n.trees=1500, interaction.depth=13, cv.folds=5, verbose=F, seed=1104) 
                                                  
                                                  
                                                  for (i in 1:ncol(X)) {
                                                    X[,i] <- as.numeric(X[,i])  # convert X to numeric values
                                                  }
                                                  
                                                  # ci model
                                                  cat("Train ci model ..."); flush.console()
                                                  mci <- CIBoost$new(maxIter=10,alpha=1,minIncr=0,gridMax=0.5,gridNum=50)
                                                  #mci$customTrain(X, clinicalSurvData, initUseAIC=TRUE)
                                                  mci$customTrain(X, clinicalSurvData)
                                                  cat("Done\n")
                                                  
                                                  # cib model
                                                  cat("Train cib model ..."); flush.console()
                                                  mcib <- ConcordanceIndexBoost$new(ntrees=100,shrinkage=1,maxdepth=6,max.stepsize=50,alpha=1)
                                                  mcib$customTrain(clinicalSurvData, X)
                                                  cat("Done\n")                 
                                                  
                                                  .self$model <- list(mcox=mcox, mgbm=mgbm, mcib=mcib, mci=mci)
                                                },
                                                
                                                
                                                
                                                customPredict = function(exprData, copyData, clinicalFeaturesData) {
                                                  rm(copyData)
                                                  
                                                  clnc <- ImputeMissingClinicalFeatures(clinicalFeaturesData, random=FALSE)                                            
                                                  clinical <- GoldiExpandClnc(clnc, use_osloval=.self$use_osloval)                                                
                                                  
                                                  meta = CreateMetageneSpace(exprs(exprData), chosenProbes = .self$chosenProbes)
                                                  ls = meta["ls",]
                                                  meta = t(apply(meta, 1, function(x){ x - median(x)}))
                                                  #===== Conditioning mesenchymal transition metagene by lymph node status
                                                  idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
                                                  mes.lymphneg = meta["mt",] * idx
                                                  mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
                                                  meta = rbind(meta, mes.lymphneg)
                                                  
                                                  idx = (clinical[,"lymph_nodes_positive"] > 3)
                                                  ls.lymphpos = ls * idx
                                                  ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
                                                  meta = rbind(meta, ls.lymphpos)
                                                  
                                                  idx = (meta["er",] < 0 & meta["erbb2", ] < 0)
                                                  ls.erneg = ls * idx
                                                  ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
                                                  meta = rbind(meta, ls.erneg)
                                                  
                                                  # prediction using secondary features different from cox clnc & expr models
                                                  X = cbind(data.frame(t(meta)), clinical)
                                                  X = X[, attr(.self$model$mcox$term, "term.labels")]
                                                  
                                                  # cox AIC model prediction                 
                                                  p1 <- predict(.self$model$mcox, X)
                                                  
                                                  # gbm model prediction                                                                                     
                                                  best.iter <- gbm.perf(.self$model$mgbm, method="cv", plot.it=FALSE)
                                                  p2 <- predict.gbm(.self$model$mgbm, X, best.iter)
                                                  
                                                  for (i in 1:ncol(X)) {
                                                    X[,i] <- as.numeric(X[,i])  # convert X to numeric values
                                                  }
                                                  
                                                  # ci model prediction
                                                  p3 <- .self$model$mci$customPredict(X)                      
                                                  
                                                  # cib model prediction
                                                  p4 <- .self$model$mcib$customPredict(X)                      
                                                  
                                                  pz <- cbind(p1,p2,p3,p4)
                                                  
                                                  .self$predictions <- pz
                                                  
                                                  pz <- scale(pz)
                                                  
                                                  if (length(.self$submit.model.index)==1) 
                                                    p <- pz[,.self$submit.model.index]
                                                  else 
                                                    p <- apply(pz[,.self$submit.model.index],1,mean)
                                                  
                                                  
                                                  return (p)
                                                }
                                                
                                              )
)


#################################################################################################


