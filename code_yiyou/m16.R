################# package #################

# install.packages("m6ALogisticModel")
# install.packages("ROCR")
# install.packages("pROC")
# install.packages("e1071")
# install.packages("doParallel")
# install.packages("party")
# install.packages("randomForest")
# install.packages("caret")
# install.packages("fitCons.UCSC.hg19")
# source("https://bioconductor.org/biocLite.R")
# biocLite("fitCons.UCSC.hg19")
# biocLite("BSgenome")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("fitCons.UCSC.hg19")
# biocLite("phastCons100way.UCSC.hg19")
# install.packages("varhandle")
# install.packages("cvAUC")
library(cvAUC)
library(varhandle)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
library(randomForest)
library(dplyr)
library(ggplot2)
library(ROCR)
library(pROC)
library(plyr)
library(dplyr)
library(e1071)
library(caret)
library(arm)
library(party)
library(data.table)
library(ggplot2)
library(m6ALogisticModel)
set.seed(998)

library(doParallel)
registerDoParallel(makePSOCKcluster(20))
setwd("~/R/x86_64-pc-linux-gnu-library/3.4")


# #To download package m6ALogisticModel, You need command: 
# devtools::install_("ZhenWei10/m6ALogisticModel")
# #some package from bioconduct,you need download them by commmand like: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome")


################## Using data in SummarizedExperiment ##############

# you need your personal path 
SE_06<-readRDS("DEseq/DEseq_06.rds")
SE_07<-readRDS("DEseq/DEseq_07.rds")
SE_07<-SE_06[which(rowRanges(SE_07)%in% rowRanges(SE_06))]
SE_08<-readRDS("DEseq/DEseq_08.rds")
SE_08<-SE_06[which(rowRanges(SE_08)%in% rowRanges(SE_06))]
SE_09<-readRDS("DEseq/DEseq_09.rds")
SE_09<-SE_06[which(rowRanges(SE_09)%in% rowRanges(SE_06))]
# get methylation sets
# each row means a methylation set



#the information of relabtive enzyme (writer, earser)
# for example, to get the information about the "METTL3-consistent"



M3_06 <- assay(SE_06[,grep("METTL3-consistent",colData(SE_06)$ID)])
M14_06 <- assay(SE_06[,grep("METTL14-consistent",colData(SE_06)$ID)])
M16_06 <- assay(SE_06[,grep("METTL16-consistent",colData(SE_06)$ID)])
fulldf_06<-NULL
fulldf_06<-data.frame(
M3<-M3_06,
M14<-M14_06,
M16<-M16_06
)
fulldf_06[is.na(fulldf_06)]<-0
fulldf_06<-fulldf_06[fulldf_06$M3==1|fulldf_06$M14==1|fulldf_06$M16==1,]
# colnames(fulldf_06)<-c("M3","M14","M16")
# write.csv(fulldf_06,"output_06.csv")
fulllist_06 <- which(M3_06==1|M14_06==1|M16_06==1)
##mettl3,14,16
# you can you use colData(SE)$ID to find which information you may need


####################### biological features generation ##################

Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

matureFE_06 <- predictors_annot(se = SE_06,
                                txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                bsgnm = Hsapiens,
                                fc = fitCons.UCSC.hg19,
                                pc = phastCons100way.UCSC.hg19,
                                struct_hybridize = Struc_hg19,
                                feature_lst = Additional_features_hg19,
                                hk_genes_list = HK_hg19_eids,
                                genes_ambiguity_method = "average")

# the each columns is a kind of bilogical features, and each row is a methylation site  
Gen_full_06 <- mcols(matureFE_06)[fulllist_06,]

Gen_full_06<-data.frame(Gen_full_06)
Gen_full_06$METTL14_TREW<-NULL
Gen_full_06$METTL16_CLIP<-NULL
Gen_full_06$METTL3_TREW<-NULL
apply(Gen_full_06, 2, function(x) length(unique(x)))
Gen_full_06$dist_sj_3_p2000<-NULL
Gen_full_06$dist_sj_5_p2000<-NULL
Gen_full_06$dist_nearest_p2000<-NULL
Gen_full_06$dist_nearest_p200<-NULL
Gen_full_06$AAACA<-NULL
Gen_full_06$AAACT<-NULL
Gen_full_06$AAACC<-NULL
Gen_full_06$AGACA<-NULL
Gen_full_06$AGACT<-NULL
Gen_full_06$AGACC<-NULL
Gen_full_06$GAACA<-NULL
Gen_full_06$GAACT<-NULL
Gen_full_06$GAACC<-NULL
Gen_full_06$GGACA<-NULL
Gen_full_06$GGACT<-NULL
Gen_full_06$GGACC<-NULL
Gen_full3_06<-Gen_full_06
Gen_full14_06<-Gen_full_06
Gen_full16_06<-Gen_full_06
Gen_full3_06$M3<-factor(fulldf_06$M3,labels=c("X0","X1"))

Gen_full14_06$M14<-factor(fulldf_06$M14,labels=c("X0","X1"))

Gen_full16_06$M16<-factor(fulldf_06$M16,labels=c("X0","X1"))


############################## sequence derived features ###############

# the sequence features from MethyRNA 

## get length of the sequences is 41 bp with the m6A motif in the center
source("Binary.R")
seq_test_06 <- as.character(DNAStringSet(Views(Hsapiens,rowRanges(SE_06)+20)))
fullseq_06<-seq_test_06[fulllist_06]
## encoding 

fullseqbin_06<-data.frame(alter_chemicalNF(fullseq_06))
apply(fullseqbin_06, 2, function(x) length(unique(x)))
fullseqbin_06$X4<-NULL
fullseqbin_06$X73<-NULL
fullseqbin_06$X77<-NULL
fullseqbin_06$X81<-NULL
fullseqbin_06$X82<-NULL
fullseqbin_06$X83<-NULL
fullseqbin_06$X85<-NULL
fullseqbin_06$X86<-NULL
fullseqbin_06$X87<-NULL
seq_full3_06<-data.frame(fullseqbin_06)
seq_full14_06<-data.frame(fullseqbin_06)
seq_full16_06<-data.frame(fullseqbin_06)
seq_full3_06$M3<-factor(fulldf_06$M3,labels = c("X0","X1"))
seq_full14_06$M14<-factor(fulldf_06$M14,labels=c("X0","X1"))
seq_full16_06$M16<-factor(fulldf_06$M16,labels=c("X0","X1"))


full_full16_06<-merge(fullseqbin_06,Gen_full16_06,by=0)
full_full16_06$Row.names<-NULL
full_inTraining16_06 <- createDataPartition(full_full16_06$M16, p = .75, list = FALSE)
full_training16_06 <- full_full16_06[ full_inTraining16_06,]
full_testing16_06  <- full_full16_06[-full_inTraining16_06,]
full_up_train16_06 <- upSample(x = full_training16_06[,-ncol(full_training16_06)],
                               y = full_training16_06$M16)    
colnames(full_up_train16_06)[colnames(full_up_train16_06)=="Class"] <- "M16"
seq_inTraining16_06 <- createDataPartition(seq_full16_06$M16, p = .75, list = FALSE)
seq_training16_06 <- seq_full16_06[ seq_inTraining16_06,]
seq_testing16_06  <- seq_full16_06[-seq_inTraining16_06,]
seq_up_train16_06 <- upSample(x = seq_training16_06[, -ncol(seq_training16_06)],
                              y = seq_training16_06$M16)    
colnames(seq_up_train16_06)[colnames(seq_up_train16_06)=="Class"] <- "M16"


Gen_inTraining16_06 <- createDataPartition(Gen_full16_06$M16, p = .75, list = FALSE)
Gen_training16_06 <- Gen_full16_06[ Gen_inTraining16_06,]
Gen_testing16_06  <- Gen_full16_06[-Gen_inTraining16_06,]
Gen_up_train16_06 <- upSample(x = Gen_training16_06[, -ncol(Gen_training16_06)],
                              y = Gen_training16_06$M16)    
colnames(Gen_up_train16_06)[colnames(Gen_up_train16_06)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_06 <- train(M16 ~ ., data = Gen_up_train16_06, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
seq_svmFit16_06 <- train(M16 ~ ., data = seq_up_train16_06, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
full_svmFit16_06 <- train(M16 ~ ., data = full_up_train16_06, 
                          method = "svmRadial", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")

Gen_RFFit3_06 <- train(M3 ~ ., data = Gen_up_train3_06, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
seq_RFFit3_06 <- train(M3 ~ ., data = seq_up_train3_06, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
full_RFFit3_06 <- train(M3 ~ ., data = full_up_train3_06, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
Gen_RFFit14_06 <- train(M14 ~ ., data = Gen_up_train14_06, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit14_06 <- train(M14 ~ ., data = seq_up_train14_06, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit14_06 <- train(M14 ~ ., data = full_up_train14_06, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_RFFit16_06 <- train(M16 ~ ., data = Gen_up_train16_06, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit16_06 <- train(M16 ~ ., data = seq_up_train16_06, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit16_06 <- train(M16 ~ ., data = full_up_train16_06, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_testing3_06$svmpred<-predict(Gen_svmFit3_06,Gen_testing3_06,type = "prob")[,2]
Gen_testing3_06$RFpred<-predict(Gen_RFFit3_06,Gen_testing3_06,type = "prob")[,2]
seq_testing3_06$svmpred<-predict(seq_svmFit3_06,seq_testing3_06,type = "prob")[,2]
seq_testing3_06$RFpred<-predict(seq_RFFit3_06,seq_testing3_06,type = "prob")[,2]
full_testing3_06$svmpred<-predict(full_svmFit3_06,full_testing3_06,type = "prob")[,2]
full_testing3_06$RFpred<-predict(full_RFFit3_06,full_testing3_06,type = "prob")[,2]

Gen_testing14_06$svmpred<-predict(Gen_svmFit14_06,Gen_testing14_06,type = "prob")[,2]
Gen_testing14_06$RFpred<-predict(Gen_RFFit14_06,Gen_testing14_06,type = "prob")[,2]
seq_testing14_06$svmpred<-predict(seq_svmFit14_06,seq_testing14_06,type = "prob")[,2]
seq_testing14_06$RFpred<-predict(seq_RFFit14_06,seq_testing14_06,type = "prob")[,2]
full_testing14_06$svmpred<-predict(full_svmFit14_06,full_testing14_06,type = "prob")[,2]
full_testing14_06$RFpred<-predict(full_RFFit14_06,full_testing14_06,type = "prob")[,2]


Gen_testing16_06$svmpred<-predict(Gen_svmFit16_06,Gen_testing16_06,type = "prob")[,2]
Gen_testing16_06$RFpred<-predict(Gen_RFFit16_06,Gen_testing16_06,type = "prob")[,2]
seq_testing16_06$svmpred<-predict(seq_svmFit16_06,seq_testing16_06,type = "prob")[,2]
seq_testing16_06$RFpred<-predict(seq_RFFit16_06,seq_testing16_06,type = "prob")[,2]
full_testing16_06$svmpred<-predict(full_svmFit16_06,full_testing16_06,type = "prob")[,2]
full_testing16_06$RFpred<-predict(full_RFFit16_06,full_testing16_06,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")




Gen_imp3_06<-varImp(Gen_RFFit3_06)$importance
Gen_imp3_06<-Gen_imp3_06[order(-Gen_imp3_06$Overall),,drop=FALSE]
Gen_imp3v_06<-as.vector(rownames(Gen_imp3_06))
seq_imp3_06<-varImp(seq_RFFit3_06)$importance
seq_imp3_06<-seq_imp3_06[order(-seq_imp3_06$Overall),,drop=FALSE]
seq_imp3v_06<-as.vector(rownames(seq_imp3_06))
full_imp3_06<-varImp(full_RFFit3_06)$importance
full_imp3_06<-full_imp3_06[order(-full_imp3_06$Overall),,drop=FALSE]
full_imp3v_06<-as.vector(rownames(full_imp3_06))

Gen_imp14_06<-varImp(Gen_RFFit14_06)$importance
Gen_imp14_06<-Gen_imp14_06[order(-Gen_imp14_06$Overall),,drop=FALSE]
Gen_imp14v_06<-as.vector(rownames(Gen_imp14_06))
seq_imp14_06<-varImp(seq_RFFit14_06)$importance
seq_imp14_06<-seq_imp14_06[order(-seq_imp14_06$Overall),,drop=FALSE]
seq_imp14v_06<-as.vector(rownames(seq_imp14_06))
full_imp14_06<-varImp(full_RFFit14_06)$importance
full_imp14_06<-full_imp14_06[order(-full_imp14_06$Overall),,drop=FALSE]
full_imp14v_06<-as.vector(rownames(full_imp14_06))

Gen_imp16_06<-varImp(Gen_RFFit16_06)$importance
Gen_imp16_06<-Gen_imp16_06[order(-Gen_imp16_06$Overall),,drop=FALSE]
Gen_imp16v_06<-as.vector(rownames(Gen_imp16_06))
seq_imp16_06<-varImp(seq_RFFit16_06)$importance
seq_imp16_06<-seq_imp16_06[order(-seq_imp16_06$Overall),,drop=FALSE]
seq_imp16v_06<-as.vector(rownames(seq_imp16_06))
full_imp16_06<-varImp(full_RFFit16_06)$importance
full_imp16_06<-full_imp16_06[order(-full_imp16_06$Overall),,drop=FALSE]
full_imp16v_06<-as.vector(rownames(full_imp16_06))

AUC_Gen_RF_16v_06<-c(0,0)
AUC_seq_RF_16v_06<-c(0,0)
AUC_full_RF_16v_06<-c(0,0)
for(i in 1:50){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_06)%in%Gen_imp16v_06[1:i]==TRUE)
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_06[,c(top_selected,ncol(Gen_up_train16_06))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_Gen_RF_16v_06[i]<-max(Gen_i16_RFfit$results$ROC)}

for(i in 1:155){
  # setcolorder(seq_up_train16, imp16v)
  top_selected <- which(colnames(seq_up_train16_06)%in%seq_imp16v_06[1:i]==TRUE)
  selected_matrix <- seq_up_train16_06[,c(top_selected,ncol(seq_up_train16_06))]
  seq_i16_RFfit<-train(M16 ~ ., data = seq_up_train16_06[,c(top_selected,ncol(seq_up_train16_06))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_seq_RF_16v_06[i]<-max(seq_i16_RFfit$results$ROC)}

for(i in 1:205){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_06)%in%full_imp16v_06[1:i]==TRUE)
  selected_matrix <- full_up_train16_06[,c(top_selected,ncol(full_up_train16_06))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_06[,c(top_selected,ncol(full_up_train16_06))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_06[i]<-max(full_i16_RFfit$results$ROC)}










top_selected <- which(colnames(full_up_train16_06)%in%full_imp16v_06[1:which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))]==TRUE)
full_16_RF_final_06<-train(M16 ~ ., data = full_up_train16_06[,c(top_selected,ncol(full_up_train16_06))], 
                           method = "rf", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
full_16_svm_final_06<-train(M16 ~ ., data = full_up_train16_06[,c(top_selected,ncol(full_up_train16_06))], 
                            method = "svmRadial", 
                            trControl = fitControl, 
                            preProc=c("center","scale"),
                            metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_06)%in%Gen_imp16v_06[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))]==TRUE)
Gen_16_RF_final_06<-train(M16 ~ ., data = Gen_up_train16_06[,c(top_selected,ncol(Gen_up_train16_06))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
Gen_16_svm_final_06<-train(M16 ~ ., data = Gen_up_train16_06[,c(top_selected,ncol(Gen_up_train16_06))], 
                           method = "svmRadialSigma", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
top_selected <- which(colnames(seq_up_train16_06)%in%seq_imp16v_06[1:which(AUC_seq_RF_16v_06==max(AUC_seq_RF_16v_06))]==TRUE)
seq_16_RF_final_06<-train(M16 ~ ., data = seq_up_train16_06[,c(top_selected,ncol(seq_up_train16_06))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
seq_16_svm_final_06<-train(M16 ~ ., data = seq_up_train16_06[,c(top_selected,ncol(seq_up_train16_06))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")


AUC_full_RF_16_final_06<-max(full_16_RF_final_06$results$ROC)
AUC_Gen_RF_16_final_06<-max(Gen_16_RF_final_06$results$ROC)
AUC_seq_RF_16_final_06<-max(seq_16_RF_final_06$results$ROC)
AUC_full_svm_16_final_06<-max(full_16_svm_final_06$results$ROC)
AUC_Gen_svm_16_final_06<-max(Gen_16_svm_final_06$results$ROC)
AUC_seq_svm_16_final_06<-max(seq_16_svm_final_06$results$ROC)




AUC_Gen_RF_16_df_06<-data.frame(AUC=AUC_Gen_RF_16v_06,number=c(1:50))
AUC_seq_RF_16_df_06<-data.frame(AUC=AUC_seq_RF_16v_06,number=c(1:155))
AUC_full_RF_16_df_06<-data.frame(AUC=AUC_full_RF_16v_06,number=c(1:205))

plot(AUC~ number,data=AUC_Gen_RF_16_df_06)
plot(AUC~ number,data=AUC_seq_RF_16_df_06)
plot(AUC~ number,data=AUC_full_RF_16_df_06)



















































































































M3_07 <- assay(SE_07[,grep("METTL3-consistent",colData(SE_07)$ID)])
M14_07 <- assay(SE_07[,grep("METTL14-consistent",colData(SE_07)$ID)])
M16_07 <- assay(SE_07[,grep("METTL16-consistent",colData(SE_07)$ID)])
fulldf_07<-data.frame(
  M3<-M3_07,
  M14<-M14_07,
  M16<-M16_07
)
fulldf_07[is.na(fulldf_07)]<-0
fulldf_07<-fulldf_07[fulldf_07$M3==1|fulldf_07$M14==1|fulldf_07$M16==1,]

fulllist_07 <- which(M3_07==1|M14_07==1|M16_07==1)
##mettl3,14,16
# you can you use colData(SE)$ID to find which information you may need


####################### biological features generation ##################

Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

matureFE_07 <- predictors_annot(se = SE_07,
                                txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                bsgnm = Hsapiens,
                                fc = fitCons.UCSC.hg19,
                                pc = phastCons100way.UCSC.hg19,
                                struct_hybridize = Struc_hg19,
                                feature_lst = Additional_features_hg19,
                                hk_genes_list = HK_hg19_eids,
                                genes_ambiguity_method = "average")

# the each columns is a kind of bilogical features, and each row is a methylation site  
Gen_full_07 <- mcols(matureFE_07)[fulllist_07,]

Gen_full_07<-data.frame(Gen_full_07)
Gen_full_07$METTL14_TREW<-NULL
Gen_full_07$METTL16_CLIP<-NULL
Gen_full_07$METTL3_TREW<-NULL
apply(Gen_full_07, 2, function(x) length(unique(x)))
Gen_full_07$dist_sj_3_p2000<-NULL
Gen_full_07$dist_sj_5_p2000<-NULL
Gen_full_07$dist_nearest_p2000<-NULL
Gen_full_07$dist_nearest_p200<-NULL
Gen_full_07$AAACA<-NULL
Gen_full_07$AAACT<-NULL
Gen_full_07$AAACC<-NULL
Gen_full_07$AGACA<-NULL
Gen_full_07$AGACT<-NULL
Gen_full_07$AGACC<-NULL
Gen_full_07$GAACA<-NULL
Gen_full_07$GAACT<-NULL
Gen_full_07$GAACC<-NULL
Gen_full_07$GGACA<-NULL
Gen_full_07$GGACT<-NULL
Gen_full_07$GGACC<-NULL
Gen_full3_07<-Gen_full_07
Gen_full14_07<-Gen_full_07
Gen_full16_07<-Gen_full_07
Gen_full3_07$M3<-factor(fulldf_07$M3,labels=c("X0","X1"))

Gen_full14_07$M14<-factor(fulldf_07$M14,labels=c("X0","X1"))

Gen_full16_07$M16<-factor(fulldf_07$M16,labels=c("X0","X1"))


############################## sequence derived features ###############

# the sequence features from MethyRNA 

## get length of the sequences is 41 bp with the m6A motif in the center
source("Binary.R")
seq_test_07 <- as.character(DNAStringSet(Views(Hsapiens,rowRanges(SE_07)+20)))
fullseq_07<-seq_test_07[fulllist_07]
## encoding 


fullseqbin_07<-data.frame(alter_chemicalNF(fullseq_07))
apply(fullseqbin_07, 2, function(x) length(unique(x)))
fullseqbin_07$X4<-NULL
fullseqbin_07$X73<-NULL
fullseqbin_07$X77<-NULL
fullseqbin_07$X81<-NULL
fullseqbin_07$X82<-NULL
fullseqbin_07$X83<-NULL
fullseqbin_07$X85<-NULL
fullseqbin_07$X86<-NULL
fullseqbin_07$X87<-NULL
seq_full3_07<-data.frame(fullseqbin_07)
seq_full14_07<-data.frame(fullseqbin_07)
seq_full16_07<-data.frame(fullseqbin_07)
seq_full3_07$M3<-factor(fulldf_07$M3,labels = c("X0","X1"))
seq_full14_07$M14<-factor(fulldf_07$M14,labels=c("X0","X1"))
seq_full16_07$M16<-factor(fulldf_07$M16,labels=c("X0","X1"))


full_full16_07<-merge(fullseqbin_07,Gen_full16_07,by=0)
full_full16_07$Row.names<-NULL
full_inTraining16_07 <- createDataPartition(full_full16_07$M16, p = .75, list = FALSE)
full_training16_07 <- full_full16_07[ full_inTraining16_07,]
full_testing16_07  <- full_full16_07[-full_inTraining16_07,]
full_up_train16_07 <- upSample(x = full_training16_07[,-ncol(full_training16_07)],
                               y = full_training16_07$M16)    
colnames(full_up_train16_07)[colnames(full_up_train16_07)=="Class"] <- "M16"
seq_inTraining16_07 <- createDataPartition(seq_full16_07$M16, p = .75, list = FALSE)
seq_training16_07 <- seq_full16_07[ seq_inTraining16_07,]
seq_testing16_07  <- seq_full16_07[-seq_inTraining16_07,]
seq_up_train16_07 <- upSample(x = seq_training16_07[, -ncol(seq_training16_07)],
                              y = seq_training16_07$M16)    
colnames(seq_up_train16_07)[colnames(seq_up_train16_07)=="Class"] <- "M16"


Gen_inTraining16_07 <- createDataPartition(Gen_full16_07$M16, p = .75, list = FALSE)
Gen_training16_07 <- Gen_full16_07[ Gen_inTraining16_07,]
Gen_testing16_07  <- Gen_full16_07[-Gen_inTraining16_07,]
Gen_up_train16_07 <- upSample(x = Gen_training16_07[, -ncol(Gen_training16_07)],
                              y = Gen_training16_07$M16)    
colnames(Gen_up_train16_07)[colnames(Gen_up_train16_07)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_07 <- train(M16 ~ ., data = Gen_up_train16_07, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
seq_svmFit16_07 <- train(M16 ~ ., data = seq_up_train16_07, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
full_svmFit16_07 <- train(M16 ~ ., data = full_up_train16_07, 
                          method = "svmRadial", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")

Gen_RFFit3_07 <- train(M3 ~ ., data = Gen_up_train3_07, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
seq_RFFit3_07 <- train(M3 ~ ., data = seq_up_train3_07, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
full_RFFit3_07 <- train(M3 ~ ., data = full_up_train3_07, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
Gen_RFFit14_07 <- train(M14 ~ ., data = Gen_up_train14_07, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit14_07 <- train(M14 ~ ., data = seq_up_train14_07, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit14_07 <- train(M14 ~ ., data = full_up_train14_07, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_RFFit16_07 <- train(M16 ~ ., data = Gen_up_train16_07, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit16_07 <- train(M16 ~ ., data = seq_up_train16_07, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit16_07 <- train(M16 ~ ., data = full_up_train16_07, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_testing3_07$svmpred<-predict(Gen_svmFit3_07,Gen_testing3_07,type = "prob")[,2]
Gen_testing3_07$RFpred<-predict(Gen_RFFit3_07,Gen_testing3_07,type = "prob")[,2]
seq_testing3_07$svmpred<-predict(seq_svmFit3_07,seq_testing3_07,type = "prob")[,2]
seq_testing3_07$RFpred<-predict(seq_RFFit3_07,seq_testing3_07,type = "prob")[,2]
full_testing3_07$svmpred<-predict(full_svmFit3_07,full_testing3_07,type = "prob")[,2]
full_testing3_07$RFpred<-predict(full_RFFit3_07,full_testing3_07,type = "prob")[,2]

Gen_testing14_07$svmpred<-predict(Gen_svmFit14_07,Gen_testing14_07,type = "prob")[,2]
Gen_testing14_07$RFpred<-predict(Gen_RFFit14_07,Gen_testing14_07,type = "prob")[,2]
seq_testing14_07$svmpred<-predict(seq_svmFit14_07,seq_testing14_07,type = "prob")[,2]
seq_testing14_07$RFpred<-predict(seq_RFFit14_07,seq_testing14_07,type = "prob")[,2]
full_testing14_07$svmpred<-predict(full_svmFit14_07,full_testing14_07,type = "prob")[,2]
full_testing14_07$RFpred<-predict(full_RFFit14_07,full_testing14_07,type = "prob")[,2]


Gen_testing16_07$svmpred<-predict(Gen_svmFit16_07,Gen_testing16_07,type = "prob")[,2]
Gen_testing16_07$RFpred<-predict(Gen_RFFit16_07,Gen_testing16_07,type = "prob")[,2]
seq_testing16_07$svmpred<-predict(seq_svmFit16_07,seq_testing16_07,type = "prob")[,2]
seq_testing16_07$RFpred<-predict(seq_RFFit16_07,seq_testing16_07,type = "prob")[,2]
full_testing16_07$svmpred<-predict(full_svmFit16_07,full_testing16_07,type = "prob")[,2]
full_testing16_07$RFpred<-predict(full_RFFit16_07,full_testing16_07,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")




Gen_imp3_07<-varImp(Gen_RFFit3_07)$importance
Gen_imp3_07<-Gen_imp3_07[order(-Gen_imp3_07$Overall),,drop=FALSE]
Gen_imp3v_07<-as.vector(rownames(Gen_imp3_07))
seq_imp3_07<-varImp(seq_RFFit3_07)$importance
seq_imp3_07<-seq_imp3_07[order(-seq_imp3_07$Overall),,drop=FALSE]
seq_imp3v_07<-as.vector(rownames(seq_imp3_07))
full_imp3_07<-varImp(full_RFFit3_07)$importance
full_imp3_07<-full_imp3_07[order(-full_imp3_07$Overall),,drop=FALSE]
full_imp3v_07<-as.vector(rownames(full_imp3_07))

Gen_imp14_07<-varImp(Gen_RFFit14_07)$importance
Gen_imp14_07<-Gen_imp14_07[order(-Gen_imp14_07$Overall),,drop=FALSE]
Gen_imp14v_07<-as.vector(rownames(Gen_imp14_07))
seq_imp14_07<-varImp(seq_RFFit14_07)$importance
seq_imp14_07<-seq_imp14_07[order(-seq_imp14_07$Overall),,drop=FALSE]
seq_imp14v_07<-as.vector(rownames(seq_imp14_07))
full_imp14_07<-varImp(full_RFFit14_07)$importance
full_imp14_07<-full_imp14_07[order(-full_imp14_07$Overall),,drop=FALSE]
full_imp14v_07<-as.vector(rownames(full_imp14_07))

Gen_imp16_07<-varImp(Gen_RFFit16_07)$importance
Gen_imp16_07<-Gen_imp16_07[order(-Gen_imp16_07$Overall),,drop=FALSE]
Gen_imp16v_07<-as.vector(rownames(Gen_imp16_07))
seq_imp16_07<-varImp(seq_RFFit16_07)$importance
seq_imp16_07<-seq_imp16_07[order(-seq_imp16_07$Overall),,drop=FALSE]
seq_imp16v_07<-as.vector(rownames(seq_imp16_07))
full_imp16_07<-varImp(full_RFFit16_07)$importance
full_imp16_07<-full_imp16_07[order(-full_imp16_07$Overall),,drop=FALSE]
full_imp16v_07<-as.vector(rownames(full_imp16_07))

AUC_Gen_RF_16v_07<-c(0,0)
AUC_seq_RF_16v_07<-c(0,0)
AUC_full_RF_16v_07<-c(0,0)

for(i in 1:50){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_07)%in%Gen_imp16v_07[1:i]==TRUE)
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_07[,c(top_selected,ncol(Gen_up_train16_07))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_Gen_RF_16v_07[i]<-max(Gen_i16_RFfit$results$ROC)}

for(i in 1:155){
  # setcolorder(seq_up_train16, imp16v)
  top_selected <- which(colnames(seq_up_train16_07)%in%seq_imp16v_07[1:i]==TRUE)
  selected_matrix <- seq_up_train16_07[,c(top_selected,ncol(seq_up_train16_07))]
  seq_i16_RFfit<-train(M16 ~ ., data = seq_up_train16_07[,c(top_selected,ncol(seq_up_train16_07))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_seq_RF_16v_07[i]<-max(seq_i16_RFfit$results$ROC)}

for(i in 1:205){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_07)%in%full_imp16v_07[1:i]==TRUE)
  selected_matrix <- full_up_train16_07[,c(top_selected,ncol(full_up_train16_07))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_07[,c(top_selected,ncol(full_up_train16_07))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_07[i]<-max(full_i16_RFfit$results$ROC)}










top_selected <- which(colnames(full_up_train16_07)%in%full_imp16v_07[1:which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))]==TRUE)
full_16_RF_final_07<-train(M16 ~ ., data = full_up_train16_07[,c(top_selected,ncol(full_up_train16_07))], 
                           method = "rf", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
full_16_svm_final_07<-train(M16 ~ ., data = full_up_train16_07[,c(top_selected,ncol(full_up_train16_07))], 
                            method = "svmRadial", 
                            trControl = fitControl, 
                            preProc=c("center","scale"),
                            metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_07)%in%Gen_imp16v_07[1:which(AUC_Gen_RF_16v_07==max(AUC_Gen_RF_16v_07))]==TRUE)
Gen_16_RF_final_07<-train(M16 ~ ., data = Gen_up_train16_07[,c(top_selected,ncol(Gen_up_train16_07))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
Gen_16_svm_final_07<-train(M16 ~ ., data = Gen_up_train16_07[,c(top_selected,ncol(Gen_up_train16_07))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
top_selected <- which(colnames(seq_up_train16_07)%in%seq_imp16v_07[1:which(AUC_seq_RF_16v_07==max(AUC_seq_RF_16v_07))]==TRUE)
seq_16_RF_final_07<-train(M16 ~ ., data = seq_up_train16_07[,c(top_selected,ncol(seq_up_train16_07))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
seq_16_svm_final_07<-train(M16 ~ ., data = seq_up_train16_07[,c(top_selected,ncol(seq_up_train16_07))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")

AUC_full_RF_16_final_07<-max(full_16_RF_final_07$results$ROC)
AUC_Gen_RF_16_final_07<-max(Gen_16_RF_final_07$results$ROC)
AUC_seq_RF_16_final_07<-max(seq_16_RF_final_07$results$ROC)
AUC_full_svm_16_final_07<-max(full_16_svm_final_07$results$ROC)
AUC_Gen_svm_16_final_07<-max(Gen_16_svm_final_07$results$ROC)
AUC_seq_svm_16_final_07<-max(seq_16_svm_final_07$results$ROC)



AUC_Gen_RF_16_df_07<-data.frame(AUC=AUC_Gen_RF_16v_07,number=c(1:50))
AUC_seq_RF_16_df_07<-data.frame(AUC=AUC_seq_RF_16v_07,number=c(1:155))
AUC_full_RF_16_df_07<-data.frame(AUC=AUC_full_RF_16v_07,number=c(1:205))

plot(AUC~ number,data=AUC_Gen_RF_16_df_07)
plot(AUC~ number,data=AUC_seq_RF_16_df_07)
plot(AUC~ number,data=AUC_full_RF_16_df_07)

















M3_08 <- assay(SE_08[,grep("METTL3-consistent",colData(SE_08)$ID)])
M14_08 <- assay(SE_08[,grep("METTL14-consistent",colData(SE_08)$ID)])
M16_08 <- assay(SE_08[,grep("METTL16-consistent",colData(SE_08)$ID)])
fulldf_08<-data.frame(
  M3<-M3_08,
  M14<-M14_08,
  M16<-M16_08
)
fulldf_08[is.na(fulldf_08)]<-0
fulldf_08<-fulldf_08[fulldf_08$M3==1|fulldf_08$M14==1|fulldf_08$M16==1,]

fulllist_08 <- which(M3_08==1|M14_08==1|M16_08==1)
##mettl3,14,16
# you can you use colData(SE)$ID to find which information you may need


####################### biological features generation ##################

Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

matureFE_08 <- predictors_annot(se = SE_08,
                                txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                bsgnm = Hsapiens,
                                fc = fitCons.UCSC.hg19,
                                pc = phastCons100way.UCSC.hg19,
                                struct_hybridize = Struc_hg19,
                                feature_lst = Additional_features_hg19,
                                hk_genes_list = HK_hg19_eids,
                                genes_ambiguity_method = "average")

# the each columns is a kind of bilogical features, and each row is a methylation site  
Gen_full_08 <- mcols(matureFE_08)[fulllist_08,]

Gen_full_08<-data.frame(Gen_full_08)
Gen_full_08$METTL14_TREW<-NULL
Gen_full_08$METTL16_CLIP<-NULL
Gen_full_08$METTL3_TREW<-NULL
apply(Gen_full_08, 2, function(x) length(unique(x)))
Gen_full_08$dist_sj_3_p2000<-NULL
Gen_full_08$dist_sj_5_p2000<-NULL
Gen_full_08$dist_nearest_p2000<-NULL
Gen_full_08$dist_nearest_p200<-NULL
Gen_full_08$AAACA<-NULL
Gen_full_08$AAACT<-NULL
Gen_full_08$AAACC<-NULL
Gen_full_08$AGACA<-NULL
Gen_full_08$AGACT<-NULL
Gen_full_08$AGACC<-NULL
Gen_full_08$GAACA<-NULL
Gen_full_08$GAACT<-NULL
Gen_full_08$GAACC<-NULL
Gen_full_08$GGACA<-NULL
Gen_full_08$GGACT<-NULL
Gen_full_08$GGACC<-NULL
Gen_full3_08<-Gen_full_08
Gen_full14_08<-Gen_full_08
Gen_full16_08<-Gen_full_08
Gen_full3_08$M3<-factor(fulldf_08$M3,labels=c("X0","X1"))

Gen_full14_08$M14<-factor(fulldf_08$M14,labels=c("X0","X1"))

Gen_full16_08$M16<-factor(fulldf_08$M16,labels=c("X0","X1"))


############################## sequence derived features ###############

# the sequence features from MethyRNA 

## get length of the sequences is 41 bp with the m6A motif in the center
source("Binary.R")
seq_test_08 <- as.character(DNAStringSet(Views(Hsapiens,rowRanges(SE_08)+20)))
fullseq_08<-seq_test_08[fulllist_08]
## encoding 


fullseqbin_08<-data.frame(alter_chemicalNF(fullseq_08))
apply(fullseqbin_08, 2, function(x) length(unique(x)))
fullseqbin_08$X4<-NULL
fullseqbin_08$X73<-NULL
fullseqbin_08$X77<-NULL
fullseqbin_08$X81<-NULL
fullseqbin_08$X82<-NULL
fullseqbin_08$X83<-NULL
fullseqbin_08$X85<-NULL
fullseqbin_08$X86<-NULL
fullseqbin_08$X87<-NULL
seq_full3_08<-data.frame(fullseqbin_08)
seq_full14_08<-data.frame(fullseqbin_08)
seq_full16_08<-data.frame(fullseqbin_08)
seq_full3_08$M3<-factor(fulldf_08$M3,labels = c("X0","X1"))
seq_full14_08$M14<-factor(fulldf_08$M14,labels=c("X0","X1"))
seq_full16_08$M16<-factor(fulldf_08$M16,labels=c("X0","X1"))


full_full16_08<-merge(fullseqbin_08,Gen_full16_08,by=0)
full_full16_08$Row.names<-NULL
full_inTraining16_08 <- createDataPartition(full_full16_08$M16, p = .75, list = FALSE)
full_training16_08 <- full_full16_08[ full_inTraining16_08,]
full_testing16_08  <- full_full16_08[-full_inTraining16_08,]
full_up_train16_08 <- upSample(x = full_training16_08[,-ncol(full_training16_08)],
                               y = full_training16_08$M16)    
colnames(full_up_train16_08)[colnames(full_up_train16_08)=="Class"] <- "M16"
seq_inTraining16_08 <- createDataPartition(seq_full16_08$M16, p = .75, list = FALSE)
seq_training16_08 <- seq_full16_08[ seq_inTraining16_08,]
seq_testing16_08  <- seq_full16_08[-seq_inTraining16_08,]
seq_up_train16_08 <- upSample(x = seq_training16_08[, -ncol(seq_training16_08)],
                              y = seq_training16_08$M16)    
colnames(seq_up_train16_08)[colnames(seq_up_train16_08)=="Class"] <- "M16"


Gen_inTraining16_08 <- createDataPartition(Gen_full16_08$M16, p = .75, list = FALSE)
Gen_training16_08 <- Gen_full16_08[ Gen_inTraining16_08,]
Gen_testing16_08  <- Gen_full16_08[-Gen_inTraining16_08,]
Gen_up_train16_08 <- upSample(x = Gen_training16_08[, -ncol(Gen_training16_08)],
                              y = Gen_training16_08$M16)    
colnames(Gen_up_train16_08)[colnames(Gen_up_train16_08)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_08 <- train(M16 ~ ., data = Gen_up_train16_08, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
seq_svmFit16_08 <- train(M16 ~ ., data = seq_up_train16_08, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
full_svmFit16_08 <- train(M16 ~ ., data = full_up_train16_08, 
                          method = "svmRadial", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")

Gen_RFFit3_08 <- train(M3 ~ ., data = Gen_up_train3_08, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
seq_RFFit3_08 <- train(M3 ~ ., data = seq_up_train3_08, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
full_RFFit3_08 <- train(M3 ~ ., data = full_up_train3_08, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
Gen_RFFit14_08 <- train(M14 ~ ., data = Gen_up_train14_08, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit14_08 <- train(M14 ~ ., data = seq_up_train14_08, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit14_08 <- train(M14 ~ ., data = full_up_train14_08, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_RFFit16_08 <- train(M16 ~ ., data = Gen_up_train16_08, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit16_08 <- train(M16 ~ ., data = seq_up_train16_08, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit16_08 <- train(M16 ~ ., data = full_up_train16_08, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_testing3_08$svmpred<-predict(Gen_svmFit3_08,Gen_testing3_08,type = "prob")[,2]
Gen_testing3_08$RFpred<-predict(Gen_RFFit3_08,Gen_testing3_08,type = "prob")[,2]
seq_testing3_08$svmpred<-predict(seq_svmFit3_08,seq_testing3_08,type = "prob")[,2]
seq_testing3_08$RFpred<-predict(seq_RFFit3_08,seq_testing3_08,type = "prob")[,2]
full_testing3_08$svmpred<-predict(full_svmFit3_08,full_testing3_08,type = "prob")[,2]
full_testing3_08$RFpred<-predict(full_RFFit3_08,full_testing3_08,type = "prob")[,2]

Gen_testing14_08$svmpred<-predict(Gen_svmFit14_08,Gen_testing14_08,type = "prob")[,2]
Gen_testing14_08$RFpred<-predict(Gen_RFFit14_08,Gen_testing14_08,type = "prob")[,2]
seq_testing14_08$svmpred<-predict(seq_svmFit14_08,seq_testing14_08,type = "prob")[,2]
seq_testing14_08$RFpred<-predict(seq_RFFit14_08,seq_testing14_08,type = "prob")[,2]
full_testing14_08$svmpred<-predict(full_svmFit14_08,full_testing14_08,type = "prob")[,2]
full_testing14_08$RFpred<-predict(full_RFFit14_08,full_testing14_08,type = "prob")[,2]


Gen_testing16_08$svmpred<-predict(Gen_svmFit16_08,Gen_testing16_08,type = "prob")[,2]
Gen_testing16_08$RFpred<-predict(Gen_RFFit16_08,Gen_testing16_08,type = "prob")[,2]
seq_testing16_08$svmpred<-predict(seq_svmFit16_08,seq_testing16_08,type = "prob")[,2]
seq_testing16_08$RFpred<-predict(seq_RFFit16_08,seq_testing16_08,type = "prob")[,2]
full_testing16_08$svmpred<-predict(full_svmFit16_08,full_testing16_08,type = "prob")[,2]
full_testing16_08$RFpred<-predict(full_RFFit16_08,full_testing16_08,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")




Gen_imp3_08<-varImp(Gen_RFFit3_08)$importance
Gen_imp3_08<-Gen_imp3_08[order(-Gen_imp3_08$Overall),,drop=FALSE]
Gen_imp3v_08<-as.vector(rownames(Gen_imp3_08))
seq_imp3_08<-varImp(seq_RFFit3_08)$importance
seq_imp3_08<-seq_imp3_08[order(-seq_imp3_08$Overall),,drop=FALSE]
seq_imp3v_08<-as.vector(rownames(seq_imp3_08))
full_imp3_08<-varImp(full_RFFit3_08)$importance
full_imp3_08<-full_imp3_08[order(-full_imp3_08$Overall),,drop=FALSE]
full_imp3v_08<-as.vector(rownames(full_imp3_08))

Gen_imp14_08<-varImp(Gen_RFFit14_08)$importance
Gen_imp14_08<-Gen_imp14_08[order(-Gen_imp14_08$Overall),,drop=FALSE]
Gen_imp14v_08<-as.vector(rownames(Gen_imp14_08))
seq_imp14_08<-varImp(seq_RFFit14_08)$importance
seq_imp14_08<-seq_imp14_08[order(-seq_imp14_08$Overall),,drop=FALSE]
seq_imp14v_08<-as.vector(rownames(seq_imp14_08))
full_imp14_08<-varImp(full_RFFit14_08)$importance
full_imp14_08<-full_imp14_08[order(-full_imp14_08$Overall),,drop=FALSE]
full_imp14v_08<-as.vector(rownames(full_imp14_08))

Gen_imp16_08<-varImp(Gen_RFFit16_08)$importance
Gen_imp16_08<-Gen_imp16_08[order(-Gen_imp16_08$Overall),,drop=FALSE]
Gen_imp16v_08<-as.vector(rownames(Gen_imp16_08))
seq_imp16_08<-varImp(seq_RFFit16_08)$importance
seq_imp16_08<-seq_imp16_08[order(-seq_imp16_08$Overall),,drop=FALSE]
seq_imp16v_08<-as.vector(rownames(seq_imp16_08))
full_imp16_08<-varImp(full_RFFit16_08)$importance
full_imp16_08<-full_imp16_08[order(-full_imp16_08$Overall),,drop=FALSE]
full_imp16v_08<-as.vector(rownames(full_imp16_08))

AUC_Gen_RF_16v_08<-c(0,0)
AUC_seq_RF_16v_08<-c(0,0)
AUC_full_RF_16v_08<-c(0,0)
for(i in 1:50){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_08)%in%Gen_imp16v_08[1:i]==TRUE)
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_08[,c(top_selected,ncol(Gen_up_train16_08))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_Gen_RF_16v_08[i]<-max(Gen_i16_RFfit$results$ROC)}

for(i in 1:155){
  # setcolorder(seq_up_train16, imp16v)
  top_selected <- which(colnames(seq_up_train16_08)%in%seq_imp16v_08[1:i]==TRUE)
  selected_matrix <- seq_up_train16_08[,c(top_selected,ncol(seq_up_train16_08))]
  seq_i16_RFfit<-train(M16 ~ ., data = seq_up_train16_08[,c(top_selected,ncol(seq_up_train16_08))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_seq_RF_16v_08[i]<-max(seq_i16_RFfit$results$ROC)}

for(i in 1:205){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_08)%in%full_imp16v_08[1:i]==TRUE)
  selected_matrix <- full_up_train16_08[,c(top_selected,ncol(full_up_train16_08))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_08[,c(top_selected,ncol(full_up_train16_08))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_08[i]<-max(full_i16_RFfit$results$ROC)}










top_selected <- which(colnames(full_up_train16_08)%in%full_imp16v_08[1:which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))]==TRUE)
full_16_RF_final_08<-train(M16 ~ ., data = full_up_train16_08[,c(top_selected,ncol(full_up_train16_08))], 
                           method = "rf", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
full_16_svm_final_08<-train(M16 ~ ., data = full_up_train16_08[,c(top_selected,ncol(full_up_train16_08))], 
                            method = "svmRadial", 
                            trControl = fitControl, 
                            preProc=c("center","scale"),
                            metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_08)%in%Gen_imp16v_08[1:which(AUC_Gen_RF_16v_08==max(AUC_Gen_RF_16v_08))]==TRUE)
Gen_16_RF_final_08<-train(M16 ~ ., data = Gen_up_train16_08[,c(top_selected,ncol(Gen_up_train16_08))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
Gen_16_svm_final_08<-train(M16 ~ ., data = Gen_up_train16_08[,c(top_selected,ncol(Gen_up_train16_08))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
top_selected <- which(colnames(seq_up_train16_08)%in%seq_imp16v_08[1:which(AUC_seq_RF_16v_08==max(AUC_seq_RF_16v_08))]==TRUE)
seq_16_RF_final_08<-train(M16 ~ ., data = seq_up_train16_08[,c(top_selected,ncol(seq_up_train16_08))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
seq_16_svm_final_08<-train(M16 ~ ., data = seq_up_train16_08[,c(top_selected,ncol(seq_up_train16_08))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")


AUC_full_RF_16_final_08<-max(full_16_RF_final_08$results$ROC)
AUC_Gen_RF_16_final_08<-max(Gen_16_RF_final_08$results$ROC)
AUC_seq_RF_16_final_08<-max(seq_16_RF_final_08$results$ROC)
AUC_full_svm_16_final_08<-max(full_16_svm_final_08$results$ROC)
AUC_Gen_svm_16_final_08<-max(Gen_16_svm_final_08$results$ROC)
AUC_seq_svm_16_final_08<-max(seq_16_svm_final_08$results$ROC)




AUC_Gen_RF_16_df_08<-data.frame(AUC=AUC_Gen_RF_16v_08,number=c(1:50))
AUC_seq_RF_16_df_08<-data.frame(AUC=AUC_seq_RF_16v_08,number=c(1:155))
AUC_full_RF_16_df_08<-data.frame(AUC=AUC_full_RF_16v_08,number=c(1:205))

plot(AUC~ number,data=AUC_Gen_RF_16_df_08)
plot(AUC~ number,data=AUC_seq_RF_16_df_08)
plot(AUC~ number,data=AUC_full_RF_16_df_08)






































M3_09 <- assay(SE_09[,grep("METTL3-consistent",colData(SE_09)$ID)])
M14_09 <- assay(SE_09[,grep("METTL14-consistent",colData(SE_09)$ID)])
M16_09 <- assay(SE_09[,grep("METTL16-consistent",colData(SE_09)$ID)])
fulldf_09<-data.frame(
  M3<-M3_09,
  M14<-M14_09,
  M16<-M16_09
)
fulldf_09[is.na(fulldf_09)]<-0
fulldf_09<-fulldf_09[fulldf_09$M3==1|fulldf_09$M14==1|fulldf_09$M16==1,]

fulllist_09 <- which(M3_09==1|M14_09==1|M16_09==1)
##mettl3,14,16
# you can you use colData(SE)$ID to find which information you may need


####################### biological features generation ##################

Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

matureFE_09 <- predictors_annot(se = SE_09,
                                txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                bsgnm = Hsapiens,
                                fc = fitCons.UCSC.hg19,
                                pc = phastCons100way.UCSC.hg19,
                                struct_hybridize = Struc_hg19,
                                feature_lst = Additional_features_hg19,
                                hk_genes_list = HK_hg19_eids,
                                genes_ambiguity_method = "average")

# the each columns is a kind of bilogical features, and each row is a methylation site  
Gen_full_09 <- mcols(matureFE_09)[fulllist_09,]

Gen_full_09<-data.frame(Gen_full_09)
Gen_full_09$METTL14_TREW<-NULL
Gen_full_09$METTL16_CLIP<-NULL
Gen_full_09$METTL3_TREW<-NULL
apply(Gen_full_09, 2, function(x) length(unique(x)))
Gen_full_09$dist_sj_3_p2000<-NULL
Gen_full_09$dist_sj_5_p2000<-NULL
Gen_full_09$dist_nearest_p2000<-NULL
Gen_full_09$dist_nearest_p200<-NULL
Gen_full_09$AAACA<-NULL
Gen_full_09$AAACT<-NULL
Gen_full_09$AAACC<-NULL
Gen_full_09$AGACA<-NULL
Gen_full_09$AGACT<-NULL
Gen_full_09$AGACC<-NULL
Gen_full_09$GAACA<-NULL
Gen_full_09$GAACT<-NULL
Gen_full_09$GAACC<-NULL
Gen_full_09$GGACA<-NULL
Gen_full_09$GGACT<-NULL
Gen_full_09$GGACC<-NULL
Gen_full3_09<-Gen_full_09
Gen_full14_09<-Gen_full_09
Gen_full16_09<-Gen_full_09
Gen_full3_09$M3<-factor(fulldf_09$M3,labels=c("X0","X1"))

Gen_full14_09$M14<-factor(fulldf_09$M14,labels=c("X0","X1"))

Gen_full16_09$M16<-factor(fulldf_09$M16,labels=c("X0","X1"))


############################## sequence derived features ###############

# the sequence features from MethyRNA 

## get length of the sequences is 41 bp with the m6A motif in the center
source("Binary.R")
seq_test_09 <- as.character(DNAStringSet(Views(Hsapiens,rowRanges(SE_09)+20)))
fullseq_09<-seq_test_09[fulllist_09]
## encoding 


fullseqbin_09<-data.frame(alter_chemicalNF(fullseq_09))
apply(fullseqbin_09, 2, function(x) length(unique(x)))
fullseqbin_09$X4<-NULL
fullseqbin_09$X73<-NULL
fullseqbin_09$X77<-NULL
fullseqbin_09$X81<-NULL
fullseqbin_09$X82<-NULL
fullseqbin_09$X83<-NULL
fullseqbin_09$X85<-NULL
fullseqbin_09$X86<-NULL
fullseqbin_09$X87<-NULL
seq_full3_09<-data.frame(fullseqbin_09)
seq_full14_09<-data.frame(fullseqbin_09)
seq_full16_09<-data.frame(fullseqbin_09)
seq_full3_09$M3<-factor(fulldf_09$M3,labels = c("X0","X1"))
seq_full14_09$M14<-factor(fulldf_09$M14,labels=c("X0","X1"))
seq_full16_09$M16<-factor(fulldf_09$M16,labels=c("X0","X1"))


full_full16_09<-merge(fullseqbin_09,Gen_full16_09,by=0)
full_full16_09$Row.names<-NULL
full_inTraining16_09 <- createDataPartition(full_full16_09$M16, p = .75, list = FALSE)
full_training16_09 <- full_full16_09[ full_inTraining16_09,]
full_testing16_09  <- full_full16_09[-full_inTraining16_09,]
full_up_train16_09 <- upSample(x = full_training16_09[,-ncol(full_training16_09)],
                               y = full_training16_09$M16)    
colnames(full_up_train16_09)[colnames(full_up_train16_09)=="Class"] <- "M16"
seq_inTraining16_09 <- createDataPartition(seq_full16_09$M16, p = .75, list = FALSE)
seq_training16_09 <- seq_full16_09[ seq_inTraining16_09,]
seq_testing16_09  <- seq_full16_09[-seq_inTraining16_09,]
seq_up_train16_09 <- upSample(x = seq_training16_09[, -ncol(seq_training16_09)],
                              y = seq_training16_09$M16)    
colnames(seq_up_train16_09)[colnames(seq_up_train16_09)=="Class"] <- "M16"


Gen_inTraining16_09 <- createDataPartition(Gen_full16_09$M16, p = .75, list = FALSE)
Gen_training16_09 <- Gen_full16_09[ Gen_inTraining16_09,]
Gen_testing16_09  <- Gen_full16_09[-Gen_inTraining16_09,]
Gen_up_train16_09 <- upSample(x = Gen_training16_09[, -ncol(Gen_training16_09)],
                              y = Gen_training16_09$M16)    
colnames(Gen_up_train16_09)[colnames(Gen_up_train16_09)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_09 <- train(M16 ~ ., data = Gen_up_train16_09, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
seq_svmFit16_09 <- train(M16 ~ ., data = seq_up_train16_09, 
                         method = "svmRadial", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
full_svmFit16_09 <- train(M16 ~ ., data = full_up_train16_09, 
                          method = "svmRadial", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")

Gen_RFFit3_09 <- train(M3 ~ ., data = Gen_up_train3_09, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
seq_RFFit3_09 <- train(M3 ~ ., data = seq_up_train3_09, 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
full_RFFit3_09 <- train(M3 ~ ., data = full_up_train3_09, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
Gen_RFFit14_09 <- train(M14 ~ ., data = Gen_up_train14_09, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit14_09 <- train(M14 ~ ., data = seq_up_train14_09, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit14_09 <- train(M14 ~ ., data = full_up_train14_09, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_RFFit16_09 <- train(M16 ~ ., data = Gen_up_train16_09, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
seq_RFFit16_09 <- train(M16 ~ ., data = seq_up_train16_09, 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
full_RFFit16_09 <- train(M16 ~ ., data = full_up_train16_09, 
                         method = "rf", 
                         trControl = fitControl, 
                         preProc=c("center","scale"),
                         metric = "ROC")
Gen_testing3_09$svmpred<-predict(Gen_svmFit3_09,Gen_testing3_09,type = "prob")[,2]
Gen_testing3_09$RFpred<-predict(Gen_RFFit3_09,Gen_testing3_09,type = "prob")[,2]
seq_testing3_09$svmpred<-predict(seq_svmFit3_09,seq_testing3_09,type = "prob")[,2]
seq_testing3_09$RFpred<-predict(seq_RFFit3_09,seq_testing3_09,type = "prob")[,2]
full_testing3_09$svmpred<-predict(full_svmFit3_09,full_testing3_09,type = "prob")[,2]
full_testing3_09$RFpred<-predict(full_RFFit3_09,full_testing3_09,type = "prob")[,2]

Gen_testing14_09$svmpred<-predict(Gen_svmFit14_09,Gen_testing14_09,type = "prob")[,2]
Gen_testing14_09$RFpred<-predict(Gen_RFFit14_09,Gen_testing14_09,type = "prob")[,2]
seq_testing14_09$svmpred<-predict(seq_svmFit14_09,seq_testing14_09,type = "prob")[,2]
seq_testing14_09$RFpred<-predict(seq_RFFit14_09,seq_testing14_09,type = "prob")[,2]
full_testing14_09$svmpred<-predict(full_svmFit14_09,full_testing14_09,type = "prob")[,2]
full_testing14_09$RFpred<-predict(full_RFFit14_09,full_testing14_09,type = "prob")[,2]


Gen_testing16_09$svmpred<-predict(Gen_svmFit16_09,Gen_testing16_09,type = "prob")[,2]
Gen_testing16_09$RFpred<-predict(Gen_RFFit16_09,Gen_testing16_09,type = "prob")[,2]
seq_testing16_09$svmpred<-predict(seq_svmFit16_09,seq_testing16_09,type = "prob")[,2]
seq_testing16_09$RFpred<-predict(seq_RFFit16_09,seq_testing16_09,type = "prob")[,2]
full_testing16_09$svmpred<-predict(full_svmFit16_09,full_testing16_09,type = "prob")[,2]
full_testing16_09$RFpred<-predict(full_RFFit16_09,full_testing16_09,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")




Gen_imp3_09<-varImp(Gen_RFFit3_09)$importance
Gen_imp3_09<-Gen_imp3_09[order(-Gen_imp3_09$Overall),,drop=FALSE]
Gen_imp3v_09<-as.vector(rownames(Gen_imp3_09))
seq_imp3_09<-varImp(seq_RFFit3_09)$importance
seq_imp3_09<-seq_imp3_09[order(-seq_imp3_09$Overall),,drop=FALSE]
seq_imp3v_09<-as.vector(rownames(seq_imp3_09))
full_imp3_09<-varImp(full_RFFit3_09)$importance
full_imp3_09<-full_imp3_09[order(-full_imp3_09$Overall),,drop=FALSE]
full_imp3v_09<-as.vector(rownames(full_imp3_09))

Gen_imp14_09<-varImp(Gen_RFFit14_09)$importance
Gen_imp14_09<-Gen_imp14_09[order(-Gen_imp14_09$Overall),,drop=FALSE]
Gen_imp14v_09<-as.vector(rownames(Gen_imp14_09))
seq_imp14_09<-varImp(seq_RFFit14_09)$importance
seq_imp14_09<-seq_imp14_09[order(-seq_imp14_09$Overall),,drop=FALSE]
seq_imp14v_09<-as.vector(rownames(seq_imp14_09))
full_imp14_09<-varImp(full_RFFit14_09)$importance
full_imp14_09<-full_imp14_09[order(-full_imp14_09$Overall),,drop=FALSE]
full_imp14v_09<-as.vector(rownames(full_imp14_09))

Gen_imp16_09<-varImp(Gen_RFFit16_09)$importance
Gen_imp16_09<-Gen_imp16_09[order(-Gen_imp16_09$Overall),,drop=FALSE]
Gen_imp16v_09<-as.vector(rownames(Gen_imp16_09))
seq_imp16_09<-varImp(seq_RFFit16_09)$importance
seq_imp16_09<-seq_imp16_09[order(-seq_imp16_09$Overall),,drop=FALSE]
seq_imp16v_09<-as.vector(rownames(seq_imp16_09))
full_imp16_09<-varImp(full_RFFit16_09)$importance
full_imp16_09<-full_imp16_09[order(-full_imp16_09$Overall),,drop=FALSE]
full_imp16v_09<-as.vector(rownames(full_imp16_09))

AUC_Gen_RF_16v_09<-c(0,0)
AUC_seq_RF_16v_09<-c(0,0)
AUC_full_RF_16v_09<-c(0,0)
for(i in 1:50){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_09)%in%Gen_imp16v_09[1:i]==TRUE)
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_09[,c(top_selected,ncol(Gen_up_train16_09))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_Gen_RF_16v_09[i]<-max(Gen_i16_RFfit$results$ROC)}

for(i in 1:155){
  # setcolorder(seq_up_train16, imp16v)
  top_selected <- which(colnames(seq_up_train16_09)%in%seq_imp16v_09[1:i]==TRUE)
  selected_matrix <- seq_up_train16_09[,c(top_selected,ncol(seq_up_train16_09))]
  seq_i16_RFfit<-train(M16 ~ ., data = seq_up_train16_09[,c(top_selected,ncol(seq_up_train16_09))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_seq_RF_16v_09[i]<-max(seq_i16_RFfit$results$ROC)}

for(i in 1:205){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_09)%in%full_imp16v_09[1:i]==TRUE)
  selected_matrix <- full_up_train16_09[,c(top_selected,ncol(full_up_train16_09))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_09[,c(top_selected,ncol(full_up_train16_09))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_09[i]<-max(full_i16_RFfit$results$ROC)}










top_selected <- which(colnames(full_up_train16_09)%in%full_imp16v_09[1:which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))]==TRUE)
full_16_RF_final_09<-train(M16 ~ ., data = full_up_train16_09[,c(top_selected,ncol(full_up_train16_09))], 
                           method = "rf", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
full_16_svm_final_09<-train(M16 ~ ., data = full_up_train16_09[,c(top_selected,ncol(full_up_train16_09))], 
                            method = "svmRadial", 
                            trControl = fitControl, 
                            preProc=c("center","scale"),
                            metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_09)%in%Gen_imp16v_09[1:which(AUC_Gen_RF_16v_09==max(AUC_Gen_RF_16v_09))]==TRUE)
Gen_16_RF_final_09<-train(M16 ~ ., data = Gen_up_train16_09[,c(top_selected,ncol(Gen_up_train16_09))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
Gen_16_svm_final_09<-train(M16 ~ ., data = Gen_up_train16_09[,c(top_selected,ncol(Gen_up_train16_09))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")
top_selected <- which(colnames(seq_up_train16_09)%in%seq_imp16v_09[1:which(AUC_seq_RF_16v_09==max(AUC_seq_RF_16v_09))]==TRUE)
seq_16_RF_final_09<-train(M16 ~ ., data = seq_up_train16_09[,c(top_selected,ncol(seq_up_train16_09))], 
                          method = "rf", 
                          trControl = fitControl, 
                          preProc=c("center","scale"),
                          metric = "ROC")
seq_16_svm_final_09<-train(M16 ~ ., data = seq_up_train16_09[,c(top_selected,ncol(seq_up_train16_09))], 
                           method = "svmRadial", 
                           trControl = fitControl, 
                           preProc=c("center","scale"),
                           metric = "ROC")


AUC_full_RF_16_final_09<-max(full_16_RF_final_09$results$ROC)
AUC_Gen_RF_16_final_09<-max(Gen_16_RF_final_09$results$ROC)
AUC_seq_RF_16_final_09<-max(seq_16_RF_final_09$results$ROC)
AUC_full_svm_16_final_09<-max(full_16_svm_final_09$results$ROC)
AUC_Gen_svm_16_final_09<-max(Gen_16_svm_final_09$results$ROC)
AUC_seq_svm_16_final_09<-max(seq_16_svm_final_09$results$ROC)




AUC_Gen_RF_16_df_09<-data.frame(AUC=AUC_Gen_RF_16v_09,number=c(1:50))
AUC_seq_RF_16_df_09<-data.frame(AUC=AUC_seq_RF_16v_09,number=c(1:155))
AUC_full_RF_16_df_09<-data.frame(AUC=AUC_full_RF_16v_09,number=c(1:205))

plot(AUC~ number,data=AUC_Gen_RF_16_df_09)
plot(AUC~ number,data=AUC_seq_RF_16_df_09)
plot(AUC~ number,data=AUC_full_RF_16_df_09)

















