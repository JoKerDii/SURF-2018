

Gen_full_06_noGC<-data.frame(Gen_full_06)
Gen_full_06_noGC$GC_cont_101bp<-NULL
Gen_full_06_noGC$GC_cont_genes<-NULL
Gen_full_06_noGC$GC_cont_101bp_abs<-NULL
Gen_full_06_noGC$METTL14_TREW<-NULL
Gen_full_06_noGC$METTL16_CLIP<-NULL
Gen_full_06_noGC$METTL3_TREW<-NULL
Gen_full_06_noGC$AAACA<-NULL
Gen_full_06_noGC$AAACT<-NULL
Gen_full_06_noGC$AAACC<-NULL
Gen_full_06_noGC$AGACA<-NULL
Gen_full_06_noGC$AGACT<-NULL
Gen_full_06_noGC$AGACC<-NULL
Gen_full_06_noGC$GAACA<-NULL
Gen_full_06_noGC$GAACT<-NULL
Gen_full_06_noGC$GAACC<-NULL
Gen_full_06_noGC$GGACA<-NULL
Gen_full_06_noGC$GGACT<-NULL
Gen_full_06_noGC$GGACC<-NULL
apply(Gen_full_06_noGC, 2, function(x) length(unique(x)))
Gen_full_06_noGC$dist_sj_3_p2000<-NULL
Gen_full_06_noGC$dist_sj_5_p2000<-NULL
Gen_full_06_noGC$dist_nearest_p2000<-NULL
Gen_full_06_noGC$dist_nearest_p200<-NULL
Gen_full3_06_noGC<-Gen_full_06_noGC
Gen_full14_06_noGC<-Gen_full_06_noGC
Gen_full16_06_noGC<-Gen_full_06_noGC
Gen_full3_06_noGC$M3<-factor(fulldf_06$M3,labels=c("X0","X1"))

Gen_full14_06_noGC$M14<-factor(fulldf_06$M14,labels=c("X0","X1"))

Gen_full16_06_noGC$M16<-factor(fulldf_06$M16,labels=c("X0","X1"))



full_full16_06_noGC<-merge(fullseqbin_06,Gen_full16_06_noGC,by=0)
full_full16_06_noGC$Row.names<-NULL
full_inTraining16_06_noGC <- createDataPartition(full_full16_06_noGC$M16, p = .75, list = FALSE)
full_training16_06_noGC <- full_full16_06_noGC[ full_inTraining16_06_noGC,]
full_testing16_06_noGC  <- full_full16_06_noGC[-full_inTraining16_06_noGC,]
full_up_train16_06_noGC <- upSample(x = full_training16_06_noGC[,-ncol(full_training16_06_noGC)],
                                    y = full_training16_06_noGC$M16)    
colnames(full_up_train16_06_noGC)[colnames(full_up_train16_06_noGC)=="Class"] <- "M16"

Gen_inTraining16_06_noGC <- createDataPartition(Gen_full16_06_noGC$M16, p = .75, list = FALSE)
Gen_training16_06_noGC <- Gen_full16_06_noGC[ Gen_inTraining16_06_noGC,]
Gen_testing16_06_noGC  <- Gen_full16_06_noGC[-Gen_inTraining16_06_noGC,]
Gen_up_train16_06_noGC <- upSample(x = Gen_training16_06_noGC[, -ncol(Gen_training16_06_noGC)],
                                   y = Gen_training16_06_noGC$M16)    
colnames(Gen_up_train16_06_noGC)[colnames(Gen_up_train16_06_noGC)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_06_noGC <- train(M16 ~ ., data = Gen_up_train16_06_noGC, 
                              method = "svmRadial", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")
seq_svmFit16_06_noGC <-seq_svmFit16_06

full_svmFit16_06_noGC <- train(M16 ~ ., data = full_up_train16_06_noGC, 
                               method = "svmRadial", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")


Gen_RFFit16_06_noGC <- train(M16 ~ ., data = Gen_up_train16_06_noGC, 
                             method = "rf", 
                             trControl = fitControl, 
                             preProc=c("center","scale"),
                             metric = "ROC")
seq_RFFit16_06_noGC <- seq_RFFit16_06
full_RFFit16_06_noGC <- train(M16 ~ ., data = full_up_train16_06_noGC, 
                              method = "rf", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")

Gen_testing16_06_noGC$svmpred<-predict(Gen_svmFit16_06_noGC,Gen_testing16_06_noGC,type = "prob")[,2]
Gen_testing16_06_noGC$RFpred<-predict(Gen_RFFit16_06_noGC,Gen_testing16_06_noGC,type = "prob")[,2]
full_testing16_06_noGC$svmpred<-predict(full_svmFit16_06_noGC,full_testing16_06_noGC,type = "prob")[,2]
full_testing16_06_noGC$RFpred<-predict(full_RFFit16_06_noGC,full_testing16_06_noGC,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")





Gen_imp16_06_noGC<-varImp(Gen_RFFit16_06_noGC)$importance
Gen_imp16_06_noGC<-Gen_imp16_06_noGC[order(-Gen_imp16_06_noGC$Overall),,drop=FALSE]
Gen_imp16v_06_noGC<-as.vector(rownames(Gen_imp16_06_noGC))
seq_imp16_06_noGC<-varImp(seq_RFFit16_06_noGC)$importance
seq_imp16_06_noGC<-seq_imp16_06_noGC[order(-seq_imp16_06_noGC$Overall),,drop=FALSE]
seq_imp16v_06_noGC<-as.vector(rownames(seq_imp16_06_noGC))
full_imp16_06_noGC<-varImp(full_RFFit16_06_noGC)$importance
full_imp16_06_noGC<-full_imp16_06_noGC[order(-full_imp16_06_noGC$Overall),,drop=FALSE]
full_imp16v_06_noGC<-as.vector(rownames(full_imp16_06_noGC))



AUC_Gen_RF_16v_06_noGC<-c(0,0)
AUC_seq_RF_16v_06_noGC<-c(0,0)
AUC_full_RF_16v_06_noGC<-c(0,0)
for(i in 1:47){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_06_noGC)%in%Gen_imp16v_06_noGC[1:i]==TRUE)
  selected_matrix <- Gen_up_train16_06_noGC[,c(top_selected,ncol(Gen_up_train16_06_noGC))]
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_06_noGC[,c(top_selected,ncol(Gen_up_train16_06_noGC))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  
  AUC_Gen_RF_16v_06_noGC[i]<-max(Gen_i16_RFfit$results$ROC)}


AUC_seq_RF_16v_06_noGC<-AUC_seq_RF_16v_06

for(i in 1:202){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_06_noGC)%in%full_imp16v_06_noGC[1:i]==TRUE)
  selected_matrix <- full_up_train16_06_noGC[,c(top_selected,ncol(full_up_train16_06_noGC))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_06_noGC[,c(top_selected,ncol(full_up_train16_06_noGC))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  
  AUC_full_RF_16v_06_noGC[i]<-max(full_i16_RFfit$results$ROC)}





top_selected <- which(colnames(full_up_train16_06_noGC)%in%full_imp16v_06_noGC[1:which(AUC_full_RF_16v_06_noGC==max(AUC_full_RF_16v_06_noGC))]==TRUE)
full_16_RF_final_06_noGC<-train(M16 ~ ., data = full_up_train16_06_noGC[,c(top_selected,ncol(full_up_train16_06_noGC))], 
                                method = "rf", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")
full_16_svm_final_06_noGC<-train(M16 ~ ., data = full_up_train16_06_noGC[,c(top_selected,ncol(full_up_train16_06_noGC))], 
                                 method = "svmRadial", 
                                 trControl = fitControl, 
                                 preProc=c("center","scale"),
                                 metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_06_noGC)%in%Gen_imp16v_06_noGC[1:which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))]==TRUE)
Gen_16_RF_final_06_noGC<-train(M16 ~ ., data = Gen_up_train16_06_noGC[,c(top_selected,ncol(Gen_up_train16_06_noGC))], 
                               method = "rf", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")
Gen_16_svm_final_06_noGC<-train(M16 ~ ., data = Gen_up_train16_06_noGC[,c(top_selected,ncol(Gen_up_train16_06_noGC))], 
                                method = "svmRadial", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")

seq_16_svm_final_06_noGC<-seq_16_svm_final_06
seq_16_RF_final_06_noGC<-seq_16_RF_final_06
AUC_full_RF_16_final_06_noGC<-max(full_16_RF_final_06_noGC$results$ROC)
AUC_Gen_RF_16_final_06_noGC<-max(Gen_16_RF_final_06_noGC$results$ROC)
AUC_seq_RF_16_final_06_noGC<-max(seq_16_RF_final_06_noGC$results$ROC)
AUC_full_svm_16_final_06_noGC<-max(full_16_svm_final_06_noGC$results$ROC)
AUC_Gen_svm_16_final_06_noGC<-max(Gen_16_svm_final_06_noGC$results$ROC)
AUC_seq_svm_16_final_06_noGC<-max(seq_16_svm_final_06_noGC$results$ROC)



AUC_Gen_RF_16_df_06_noGC<-data.frame(AUC=AUC_Gen_RF_16v_06_noGC,number=c(1:47))
AUC_seq_RF_16_df_06_noGC<-data.frame(AUC=AUC_seq_RF_16v_06_noGC,number=c(1:155))
AUC_full_RF_16_df_06_noGC<-data.frame(AUC=AUC_full_RF_16v_06_noGC,number=c(1:202))

plot(AUC~ number,data=AUC_Gen_RF_16_df_06_noGC)
plot(AUC~ number,data=AUC_seq_RF_16_df_06_noGC)
plot(AUC~ number,data=AUC_full_RF_16_df_06_noGC)

































Gen_full_07_noGC<-data.frame(Gen_full_07)
Gen_full_07_noGC$GC_cont_101bp<-NULL
Gen_full_07_noGC$GC_cont_genes<-NULL
Gen_full_07_noGC$GC_cont_101bp_abs<-NULL
Gen_full_07_noGC$METTL14_TREW<-NULL
Gen_full_07_noGC$METTL16_CLIP<-NULL
Gen_full_07_noGC$METTL3_TREW<-NULL
Gen_full_07_noGC$AAACA<-NULL
Gen_full_07_noGC$AAACT<-NULL
Gen_full_07_noGC$AAACC<-NULL
Gen_full_07_noGC$AGACA<-NULL
Gen_full_07_noGC$AGACT<-NULL
Gen_full_07_noGC$AGACC<-NULL
Gen_full_07_noGC$GAACA<-NULL
Gen_full_07_noGC$GAACT<-NULL
Gen_full_07_noGC$GAACC<-NULL
Gen_full_07_noGC$GGACA<-NULL
Gen_full_07_noGC$GGACT<-NULL
Gen_full_07_noGC$GGACC<-NULL
apply(Gen_full_07_noGC, 2, function(x) length(unique(x)))
Gen_full_07_noGC$dist_sj_3_p2000<-NULL
Gen_full_07_noGC$dist_sj_5_p2000<-NULL
Gen_full_07_noGC$dist_nearest_p2000<-NULL
Gen_full_07_noGC$dist_nearest_p200<-NULL
Gen_full3_07_noGC<-Gen_full_07_noGC
Gen_full14_07_noGC<-Gen_full_07_noGC
Gen_full16_07_noGC<-Gen_full_07_noGC
Gen_full3_07_noGC$M3<-factor(fulldf_07$M3,labels=c("X0","X1"))

Gen_full14_07_noGC$M14<-factor(fulldf_07$M14,labels=c("X0","X1"))

Gen_full16_07_noGC$M16<-factor(fulldf_07$M16,labels=c("X0","X1"))



full_full16_07_noGC<-merge(fullseqbin_07,Gen_full16_07_noGC,by=0)
full_full16_07_noGC$Row.names<-NULL
full_inTraining16_07_noGC <- createDataPartition(full_full16_07_noGC$M16, p = .75, list = FALSE)
full_training16_07_noGC <- full_full16_07_noGC[ full_inTraining16_07_noGC,]
full_testing16_07_noGC  <- full_full16_07_noGC[-full_inTraining16_07_noGC,]
full_up_train16_07_noGC <- upSample(x = full_training16_07_noGC[,-ncol(full_training16_07_noGC)],
                                    y = full_training16_07_noGC$M16)    
colnames(full_up_train16_07_noGC)[colnames(full_up_train16_07_noGC)=="Class"] <- "M16"

Gen_inTraining16_07_noGC <- createDataPartition(Gen_full16_07_noGC$M16, p = .75, list = FALSE)
Gen_training16_07_noGC <- Gen_full16_07_noGC[ Gen_inTraining16_07_noGC,]
Gen_testing16_07_noGC  <- Gen_full16_07_noGC[-Gen_inTraining16_07_noGC,]
Gen_up_train16_07_noGC <- upSample(x = Gen_training16_07_noGC[, -ncol(Gen_training16_07_noGC)],
                                   y = Gen_training16_07_noGC$M16)    
colnames(Gen_up_train16_07_noGC)[colnames(Gen_up_train16_07_noGC)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_07_noGC <- train(M16 ~ ., data = Gen_up_train16_07_noGC, 
                              method = "svmRadial", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")
seq_svmFit16_07_noGC <-seq_svmFit16_07

full_svmFit16_07_noGC <- train(M16 ~ ., data = full_up_train16_07_noGC, 
                               method = "svmRadial", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")


Gen_RFFit16_07_noGC <- train(M16 ~ ., data = Gen_up_train16_07_noGC, 
                             method = "rf", 
                             trControl = fitControl, 
                             preProc=c("center","scale"),
                             metric = "ROC")
seq_RFFit16_07_noGC <- seq_RFFit16_07
full_RFFit16_07_noGC <- train(M16 ~ ., data = full_up_train16_07_noGC, 
                              method = "rf", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")

Gen_testing16_07_noGC$svmpred<-predict(Gen_svmFit16_07_noGC,Gen_testing16_07_noGC,type = "prob")[,2]
Gen_testing16_07_noGC$RFpred<-predict(Gen_RFFit16_07_noGC,Gen_testing16_07_noGC,type = "prob")[,2]
full_testing16_07_noGC$svmpred<-predict(full_svmFit16_07_noGC,full_testing16_07_noGC,type = "prob")[,2]
full_testing16_07_noGC$RFpred<-predict(full_RFFit16_07_noGC,full_testing16_07_noGC,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")





Gen_imp16_07_noGC<-varImp(Gen_RFFit16_07_noGC)$importance
Gen_imp16_07_noGC<-Gen_imp16_07_noGC[order(-Gen_imp16_07_noGC$Overall),,drop=FALSE]
Gen_imp16v_07_noGC<-as.vector(rownames(Gen_imp16_07_noGC))
seq_imp16_07_noGC<-varImp(seq_RFFit16_07_noGC)$importance
seq_imp16_07_noGC<-seq_imp16_07_noGC[order(-seq_imp16_07_noGC$Overall),,drop=FALSE]
seq_imp16v_07_noGC<-as.vector(rownames(seq_imp16_07_noGC))
full_imp16_07_noGC<-varImp(full_RFFit16_07_noGC)$importance
full_imp16_07_noGC<-full_imp16_07_noGC[order(-full_imp16_07_noGC$Overall),,drop=FALSE]
full_imp16v_07_noGC<-as.vector(rownames(full_imp16_07_noGC))



AUC_Gen_RF_16v_07_noGC<-c(0,0)
AUC_seq_RF_16v_07_noGC<-c(0,0)
AUC_full_RF_16v_07_noGC<-c(0,0)
for(i in 1:47){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_07_noGC)%in%Gen_imp16v_07_noGC[1:i]==TRUE)
  selected_matrix <- Gen_up_train16_07_noGC[,c(top_selected,ncol(Gen_up_train16_07_noGC))]
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_07_noGC[,c(top_selected,ncol(Gen_up_train16_07_noGC))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  
  AUC_Gen_RF_16v_07_noGC[i]<-max(Gen_i16_RFfit$results$ROC)}


AUC_seq_RF_16v_07_noGC<-AUC_seq_RF_16v_07

for(i in 1:202){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_07_noGC)%in%full_imp16v_07_noGC[1:i]==TRUE)
  selected_matrix <- full_up_train16_07_noGC[,c(top_selected,ncol(full_up_train16_07_noGC))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_07_noGC[,c(top_selected,ncol(full_up_train16_07_noGC))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_07_noGC[i]<-max(full_i16_RFfit$results$ROC)}





top_selected <- which(colnames(full_up_train16_07_noGC)%in%full_imp16v_07_noGC[1:which(AUC_full_RF_16v_07_noGC==max(AUC_full_RF_16v_07_noGC))]==TRUE)
full_16_RF_final_07_noGC<-train(M16 ~ ., data = full_up_train16_07_noGC[,c(top_selected,ncol(full_up_train16_07_noGC))], 
                                method = "rf", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")
full_16_svm_final_07_noGC<-train(M16 ~ ., data = full_up_train16_07_noGC[,c(top_selected,ncol(full_up_train16_07_noGC))], 
                                 method = "svmRadial", 
                                 trControl = fitControl, 
                                 preProc=c("center","scale"),
                                 metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_07_noGC)%in%Gen_imp16v_07_noGC[1:which(AUC_Gen_RF_16v_07_noGC==max(AUC_Gen_RF_16v_07_noGC))]==TRUE)
Gen_16_RF_final_07_noGC<-train(M16 ~ ., data = Gen_up_train16_07_noGC[,c(top_selected,ncol(Gen_up_train16_07_noGC))], 
                               method = "rf", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")
Gen_16_svm_final_07_noGC<-train(M16 ~ ., data = Gen_up_train16_07_noGC[,c(top_selected,ncol(Gen_up_train16_07_noGC))], 
                                method = "svmRadial", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")
seq_16_RF_final_07_noGC<-seq_16_RF_final_07
seq_16_svm_final_07_noGC<-seq_16_svm_final_07_noGC
AUC_full_RF_16_final_07_noGC<-max(full_16_RF_final_07_noGC$results$ROC)
AUC_Gen_RF_16_final_07_noGC<-max(Gen_16_RF_final_07_noGC$results$ROC)
AUC_seq_RF_16_final_07_noGC<-max(seq_16_RF_final_07_noGC$results$ROC)
AUC_full_svm_16_final_07_noGC<-max(full_16_svm_final_07_noGC$results$ROC)
AUC_Gen_svm_16_final_07_noGC<-max(Gen_16_svm_final_07_noGC$results$ROC)
AUC_seq_svm_16_final_07_noGC<-max(seq_16_svm_final_07_noGC$results$ROC)



AUC_Gen_RF_16_df_07_noGC<-data.frame(AUC=AUC_Gen_RF_16v_07_noGC,number=c(1:47))
AUC_seq_RF_16_df_07_noGC<-data.frame(AUC=AUC_seq_RF_16v_07_noGC,number=c(1:155))
AUC_full_RF_16_df_07_noGC<-data.frame(AUC=AUC_full_RF_16v_07_noGC,number=c(1:202))

plot(AUC~ number,data=AUC_Gen_RF_16_df_07_noGC)
plot(AUC~ number,data=AUC_seq_RF_16_df_07_noGC)
plot(AUC~ number,data=AUC_full_RF_16_df_07_noGC)
































Gen_full_08_noGC<-data.frame(Gen_full_08)
Gen_full_08_noGC$GC_cont_101bp<-NULL
Gen_full_08_noGC$GC_cont_genes<-NULL
Gen_full_08_noGC$GC_cont_101bp_abs<-NULL
Gen_full_08_noGC$METTL14_TREW<-NULL
Gen_full_08_noGC$METTL16_CLIP<-NULL
Gen_full_08_noGC$METTL3_TREW<-NULL
Gen_full_08_noGC$AAACA<-NULL
Gen_full_08_noGC$AAACT<-NULL
Gen_full_08_noGC$AAACC<-NULL
Gen_full_08_noGC$AGACA<-NULL
Gen_full_08_noGC$AGACT<-NULL
Gen_full_08_noGC$AGACC<-NULL
Gen_full_08_noGC$GAACA<-NULL
Gen_full_08_noGC$GAACT<-NULL
Gen_full_08_noGC$GAACC<-NULL
Gen_full_08_noGC$GGACA<-NULL
Gen_full_08_noGC$GGACT<-NULL
Gen_full_08_noGC$GGACC<-NULL
apply(Gen_full_08_noGC, 2, function(x) length(unique(x)))
Gen_full_08_noGC$dist_sj_3_p2000<-NULL
Gen_full_08_noGC$dist_sj_5_p2000<-NULL
Gen_full_08_noGC$dist_nearest_p2000<-NULL
Gen_full_08_noGC$dist_nearest_p200<-NULL
Gen_full3_08_noGC<-Gen_full_08_noGC
Gen_full14_08_noGC<-Gen_full_08_noGC
Gen_full16_08_noGC<-Gen_full_08_noGC
Gen_full3_08_noGC$M3<-factor(fulldf_08$M3,labels=c("X0","X1"))

Gen_full14_08_noGC$M14<-factor(fulldf_08$M14,labels=c("X0","X1"))

Gen_full16_08_noGC$M16<-factor(fulldf_08$M16,labels=c("X0","X1"))



full_full16_08_noGC<-merge(fullseqbin_08,Gen_full16_08_noGC,by=0)
full_full16_08_noGC$Row.names<-NULL
full_inTraining16_08_noGC <- createDataPartition(full_full16_08_noGC$M16, p = .75, list = FALSE)
full_training16_08_noGC <- full_full16_08_noGC[ full_inTraining16_08_noGC,]
full_testing16_08_noGC  <- full_full16_08_noGC[-full_inTraining16_08_noGC,]
full_up_train16_08_noGC <- upSample(x = full_training16_08_noGC[,-ncol(full_training16_08_noGC)],
                                    y = full_training16_08_noGC$M16)    
colnames(full_up_train16_08_noGC)[colnames(full_up_train16_08_noGC)=="Class"] <- "M16"

Gen_inTraining16_08_noGC <- createDataPartition(Gen_full16_08_noGC$M16, p = .75, list = FALSE)
Gen_training16_08_noGC <- Gen_full16_08_noGC[ Gen_inTraining16_08_noGC,]
Gen_testing16_08_noGC  <- Gen_full16_08_noGC[-Gen_inTraining16_08_noGC,]
Gen_up_train16_08_noGC <- upSample(x = Gen_training16_08_noGC[, -ncol(Gen_training16_08_noGC)],
                                   y = Gen_training16_08_noGC$M16)    
colnames(Gen_up_train16_08_noGC)[colnames(Gen_up_train16_08_noGC)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_08_noGC <- train(M16 ~ ., data = Gen_up_train16_08_noGC, 
                              method = "svmRadial", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")
seq_svmFit16_08_noGC <-seq_svmFit16_08

full_svmFit16_08_noGC <- train(M16 ~ ., data = full_up_train16_08_noGC, 
                               method = "svmRadial", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")


Gen_RFFit16_08_noGC <- train(M16 ~ ., data = Gen_up_train16_08_noGC, 
                             method = "rf", 
                             trControl = fitControl, 
                             preProc=c("center","scale"),
                             metric = "ROC")
seq_RFFit16_08_noGC <- seq_RFFit16_08
full_RFFit16_08_noGC <- train(M16 ~ ., data = full_up_train16_08_noGC, 
                              method = "rf", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")

Gen_testing16_08_noGC$svmpred<-predict(Gen_svmFit16_08_noGC,Gen_testing16_08_noGC,type = "prob")[,2]
Gen_testing16_08_noGC$RFpred<-predict(Gen_RFFit16_08_noGC,Gen_testing16_08_noGC,type = "prob")[,2]
full_testing16_08_noGC$svmpred<-predict(full_svmFit16_08_noGC,full_testing16_08_noGC,type = "prob")[,2]
full_testing16_08_noGC$RFpred<-predict(full_RFFit16_08_noGC,full_testing16_08_noGC,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")





Gen_imp16_08_noGC<-varImp(Gen_RFFit16_08_noGC)$importance
Gen_imp16_08_noGC<-Gen_imp16_08_noGC[order(-Gen_imp16_08_noGC$Overall),,drop=FALSE]
Gen_imp16v_08_noGC<-as.vector(rownames(Gen_imp16_08_noGC))
seq_imp16_08_noGC<-varImp(seq_RFFit16_08_noGC)$importance
seq_imp16_08_noGC<-seq_imp16_08_noGC[order(-seq_imp16_08_noGC$Overall),,drop=FALSE]
seq_imp16v_08_noGC<-as.vector(rownames(seq_imp16_08_noGC))
full_imp16_08_noGC<-varImp(full_RFFit16_08_noGC)$importance
full_imp16_08_noGC<-full_imp16_08_noGC[order(-full_imp16_08_noGC$Overall),,drop=FALSE]
full_imp16v_08_noGC<-as.vector(rownames(full_imp16_08_noGC))



AUC_Gen_RF_16v_08_noGC<-c(0,0)
AUC_seq_RF_16v_08_noGC<-c(0,0)
AUC_full_RF_16v_08_noGC<-c(0,0)
for(i in 1:47){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_08_noGC)%in%Gen_imp16v_08_noGC[1:i]==TRUE)
  selected_matrix <- Gen_up_train16_08_noGC[,c(top_selected,ncol(Gen_up_train16_08_noGC))]
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_08_noGC[,c(top_selected,ncol(Gen_up_train16_08_noGC))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_Gen_RF_16v_08_noGC[i]<-max(Gen_i16_RFfit$results$ROC)}


AUC_seq_RF_16v_08_noGC<-AUC_seq_RF_16v_08

for(i in 1:202){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_08_noGC)%in%full_imp16v_08_noGC[1:i]==TRUE)
  selected_matrix <- full_up_train16_08_noGC[,c(top_selected,ncol(full_up_train16_08_noGC))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_08_noGC[,c(top_selected,ncol(full_up_train16_08_noGC))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_08_noGC[i]<-max(full_i16_RFfit$results$ROC)}





top_selected <- which(colnames(full_up_train16_08_noGC)%in%full_imp16v_08_noGC[1:which(AUC_full_RF_16v_08_noGC==max(AUC_full_RF_16v_08_noGC))]==TRUE)
full_16_RF_final_08_noGC<-train(M16 ~ ., data = full_up_train16_08_noGC[,c(top_selected,ncol(full_up_train16_08_noGC))], 
                                method = "rf", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")
full_16_svm_final_08_noGC<-train(M16 ~ ., data = full_up_train16_08_noGC[,c(top_selected,ncol(full_up_train16_08_noGC))], 
                                 method = "svmRadial", 
                                 trControl = fitControl, 
                                 preProc=c("center","scale"),
                                 metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_08_noGC)%in%Gen_imp16v_08_noGC[1:which(AUC_Gen_RF_16v_08_noGC==max(AUC_Gen_RF_16v_08_noGC))]==TRUE)
Gen_16_RF_final_08_noGC<-train(M16 ~ ., data = Gen_up_train16_08_noGC[,c(top_selected,ncol(Gen_up_train16_08_noGC))], 
                               method = "rf", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")
Gen_16_svm_final_08_noGC<-train(M16 ~ ., data = Gen_up_train16_08_noGC[,c(top_selected,ncol(Gen_up_train16_08_noGC))], 
                                method = "svmRadial", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")

seq_16_RF_final_08_noGC <- seq_16_RF_final_08
seq_16_svm_final_08_noGC <- seq_16_svm_final_08
AUC_full_RF_16_final_08_noGC<-max(full_16_RF_final_08_noGC$results$ROC)
AUC_Gen_RF_16_final_08_noGC<-max(Gen_16_RF_final_08_noGC$results$ROC)
AUC_seq_RF_16_final_08_noGC<-max(seq_16_RF_final_08_noGC$results$ROC)
AUC_full_svm_16_final_08_noGC<-max(full_16_svm_final_08_noGC$results$ROC)
AUC_Gen_svm_16_final_08_noGC<-max(Gen_16_svm_final_08_noGC$results$ROC)
AUC_seq_svm_16_final_08_noGC<-max(seq_16_svm_final_08_noGC$results$ROC)


AUC_Gen_RF_16_df_08_noGC<-data.frame(AUC=AUC_Gen_RF_16v_08_noGC,number=c(1:47))
AUC_seq_RF_16_df_08_noGC<-data.frame(AUC=AUC_seq_RF_16v_08_noGC,number=c(1:155))
AUC_full_RF_16_df_08_noGC<-data.frame(AUC=AUC_full_RF_16v_08_noGC,number=c(1:202))

plot(AUC~ number,data=AUC_Gen_RF_16_df_08_noGC)
plot(AUC~ number,data=AUC_seq_RF_16_df_08_noGC)
plot(AUC~ number,data=AUC_full_RF_16_df_08_noGC)

































Gen_full_09_noGC<-data.frame(Gen_full_09)
Gen_full_09_noGC$GC_cont_101bp<-NULL
Gen_full_09_noGC$GC_cont_genes<-NULL
Gen_full_09_noGC$GC_cont_101bp_abs<-NULL
Gen_full_09_noGC$METTL14_TREW<-NULL
Gen_full_09_noGC$METTL16_CLIP<-NULL
Gen_full_09_noGC$METTL3_TREW<-NULL
Gen_full_09_noGC$AAACA<-NULL
Gen_full_09_noGC$AAACT<-NULL
Gen_full_09_noGC$AAACC<-NULL
Gen_full_09_noGC$AGACA<-NULL
Gen_full_09_noGC$AGACT<-NULL
Gen_full_09_noGC$AGACC<-NULL
Gen_full_09_noGC$GAACA<-NULL
Gen_full_09_noGC$GAACT<-NULL
Gen_full_09_noGC$GAACC<-NULL
Gen_full_09_noGC$GGACA<-NULL
Gen_full_09_noGC$GGACT<-NULL
Gen_full_09_noGC$GGACC<-NULL
apply(Gen_full_09_noGC, 2, function(x) length(unique(x)))
Gen_full_09_noGC$dist_sj_3_p2000<-NULL
Gen_full_09_noGC$dist_sj_5_p2000<-NULL
Gen_full_09_noGC$dist_nearest_p2000<-NULL
Gen_full_09_noGC$dist_nearest_p200<-NULL
Gen_full3_09_noGC<-Gen_full_09_noGC
Gen_full14_09_noGC<-Gen_full_09_noGC
Gen_full16_09_noGC<-Gen_full_09_noGC
Gen_full3_09_noGC$M3<-factor(fulldf_09$M3,labels=c("X0","X1"))

Gen_full14_09_noGC$M14<-factor(fulldf_09$M14,labels=c("X0","X1"))

Gen_full16_09_noGC$M16<-factor(fulldf_09$M16,labels=c("X0","X1"))



full_full16_09_noGC<-merge(fullseqbin_09,Gen_full16_09_noGC,by=0)
full_full16_09_noGC$Row.names<-NULL
full_inTraining16_09_noGC <- createDataPartition(full_full16_09_noGC$M16, p = .75, list = FALSE)
full_training16_09_noGC <- full_full16_09_noGC[ full_inTraining16_09_noGC,]
full_testing16_09_noGC  <- full_full16_09_noGC[-full_inTraining16_09_noGC,]
full_up_train16_09_noGC <- upSample(x = full_training16_09_noGC[,-ncol(full_training16_09_noGC)],
                                    y = full_training16_09_noGC$M16)    
colnames(full_up_train16_09_noGC)[colnames(full_up_train16_09_noGC)=="Class"] <- "M16"

Gen_inTraining16_09_noGC <- createDataPartition(Gen_full16_09_noGC$M16, p = .75, list = FALSE)
Gen_training16_09_noGC <- Gen_full16_09_noGC[ Gen_inTraining16_09_noGC,]
Gen_testing16_09_noGC  <- Gen_full16_09_noGC[-Gen_inTraining16_09_noGC,]
Gen_up_train16_09_noGC <- upSample(x = Gen_training16_09_noGC[, -ncol(Gen_training16_09_noGC)],
                                   y = Gen_training16_09_noGC$M16)    
colnames(Gen_up_train16_09_noGC)[colnames(Gen_up_train16_09_noGC)=="Class"] <- "M16"

fitControl <- trainControl(method = "cv",
                           number = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

Gen_svmFit16_09_noGC <- train(M16 ~ ., data = Gen_up_train16_09_noGC, 
                              method = "svmRadial", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")
seq_svmFit16_09_noGC <-seq_svmFit16_09

full_svmFit16_09_noGC <- train(M16 ~ ., data = full_up_train16_09_noGC, 
                               method = "svmRadial", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")


Gen_RFFit16_09_noGC <- train(M16 ~ ., data = Gen_up_train16_09_noGC, 
                             method = "rf", 
                             trControl = fitControl, 
                             preProc=c("center","scale"),
                             metric = "ROC")
seq_RFFit16_09_noGC <- seq_RFFit16_09
full_RFFit16_09_noGC <- train(M16 ~ ., data = full_up_train16_09_noGC, 
                              method = "rf", 
                              trControl = fitControl, 
                              preProc=c("center","scale"),
                              metric = "ROC")

Gen_testing16_09_noGC$svmpred<-predict(Gen_svmFit16_09_noGC,Gen_testing16_09_noGC,type = "prob")[,2]
Gen_testing16_09_noGC$RFpred<-predict(Gen_RFFit16_09_noGC,Gen_testing16_09_noGC,type = "prob")[,2]
full_testing16_09_noGC$svmpred<-predict(full_svmFit16_09_noGC,full_testing16_09_noGC,type = "prob")[,2]
full_testing16_09_noGC$RFpred<-predict(full_RFFit16_09_noGC,full_testing16_09_noGC,type = "prob")[,2]
#pred_Gen_svm_3<-prediction(Gen_testing3$svmpred,Gen_testing3$M3)
#ROC_Gen_svm_3<-performance(pred_Gen_svm_3,measure = "auc")





Gen_imp16_09_noGC<-varImp(Gen_RFFit16_09_noGC)$importance
Gen_imp16_09_noGC<-Gen_imp16_09_noGC[order(-Gen_imp16_09_noGC$Overall),,drop=FALSE]
Gen_imp16v_09_noGC<-as.vector(rownames(Gen_imp16_09_noGC))
seq_imp16_09_noGC<-varImp(seq_RFFit16_09_noGC)$importance
seq_imp16_09_noGC<-seq_imp16_09_noGC[order(-seq_imp16_09_noGC$Overall),,drop=FALSE]
seq_imp16v_09_noGC<-as.vector(rownames(seq_imp16_09_noGC))
full_imp16_09_noGC<-varImp(full_RFFit16_09_noGC)$importance
full_imp16_09_noGC<-full_imp16_09_noGC[order(-full_imp16_09_noGC$Overall),,drop=FALSE]
full_imp16v_09_noGC<-as.vector(rownames(full_imp16_09_noGC))



AUC_Gen_RF_16v_09_noGC<-c(0,0)
AUC_seq_RF_16v_09_noGC<-c(0,0)
AUC_full_RF_16v_09_noGC<-c(0,0)
for(i in 1:47){
  # setcolorder(Gen_up_train16, imp16v)
  top_selected <- which(colnames(Gen_up_train16_09_noGC)%in%Gen_imp16v_09_noGC[1:i]==TRUE)
  selected_matrix <- Gen_up_train16_09_noGC[,c(top_selected,ncol(Gen_up_train16_09_noGC))]
  Gen_i16_RFfit<-train(M16 ~ ., data = Gen_up_train16_09_noGC[,c(top_selected,ncol(Gen_up_train16_09_noGC))], 
                       method = "rf", 
                       trControl = fitControl, 
                       preProc=c("center","scale"),
                       metric = "ROC")
  AUC_Gen_RF_16v_09_noGC[i]<-max(Gen_i16_RFfit$results$ROC)}


AUC_seq_RF_16v_09_noGC<-AUC_seq_RF_16v_09

for(i in 1:202){
  # setcolorder(full_up_train16, imp16v)
  top_selected <- which(colnames(full_up_train16_09_noGC)%in%full_imp16v_09_noGC[1:i]==TRUE)
  selected_matrix <- full_up_train16_09_noGC[,c(top_selected,ncol(full_up_train16_09_noGC))]
  full_i16_RFfit<-train(M16 ~ ., data = full_up_train16_09_noGC[,c(top_selected,ncol(full_up_train16_09_noGC))], 
                        method = "rf", 
                        trControl = fitControl, 
                        preProc=c("center","scale"),
                        metric = "ROC")
  AUC_full_RF_16v_09_noGC[i]<-max(full_i16_RFfit$results$ROC)}





top_selected <- which(colnames(full_up_train16_09_noGC)%in%full_imp16v_09_noGC[1:which(AUC_full_RF_16v_09_noGC==max(AUC_full_RF_16v_09_noGC))]==TRUE)
full_16_RF_final_09_noGC<-train(M16 ~ ., data = full_up_train16_09_noGC[,c(top_selected,ncol(full_up_train16_09_noGC))], 
                                method = "rf", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")
full_16_svm_final_09_noGC<-train(M16 ~ ., data = full_up_train16_09_noGC[,c(top_selected,ncol(full_up_train16_09_noGC))], 
                                 method = "svmRadial", 
                                 trControl = fitControl, 
                                 preProc=c("center","scale"),
                                 metric = "ROC")
top_selected <- which(colnames(Gen_up_train16_09_noGC)%in%Gen_imp16v_09_noGC[1:which(AUC_Gen_RF_16v_09_noGC==max(AUC_Gen_RF_16v_09_noGC))]==TRUE)
Gen_16_RF_final_09_noGC<-train(M16 ~ ., data = Gen_up_train16_09_noGC[,c(top_selected,ncol(Gen_up_train16_09_noGC))], 
                               method = "rf", 
                               trControl = fitControl, 
                               preProc=c("center","scale"),
                               metric = "ROC")
Gen_16_svm_final_09_noGC<-train(M16 ~ ., data = Gen_up_train16_09_noGC[,c(top_selected,ncol(Gen_up_train16_09_noGC))], 
                                method = "svmRadial", 
                                trControl = fitControl, 
                                preProc=c("center","scale"),
                                metric = "ROC")

seq_16_RF_final_09_noGC<-seq_16_RF_final_09
seq_16_svm_final_09_noGC<-seq_16_svm_final_09
AUC_full_RF_16_final_09_noGC<-max(full_16_RF_final_09_noGC$results$ROC)
AUC_Gen_RF_16_final_09_noGC<-max(Gen_16_RF_final_09_noGC$results$ROC)
AUC_seq_RF_16_final_09_noGC<-max(seq_16_RF_final_09_noGC$results$ROC)
AUC_full_svm_16_final_09_noGC<-max(full_16_svm_final_09_noGC$results$ROC)
AUC_Gen_svm_16_final_09_noGC<-max(Gen_16_svm_final_09_noGC$results$ROC)
AUC_seq_svm_16_final_09_noGC<-max(seq_16_svm_final_09_noGC$results$ROC)


AUC_Gen_RF_16_df_09_noGC<-data.frame(AUC=AUC_Gen_RF_16v_09_noGC,number=c(1:47))
AUC_seq_RF_16_df_09_noGC<-data.frame(AUC=AUC_seq_RF_16v_09_noGC,number=c(1:155))
AUC_full_RF_16_df_09_noGC<-data.frame(AUC=AUC_full_RF_16v_09_noGC,number=c(1:202))

plot(AUC~ number,data=AUC_Gen_RF_16_df_09_noGC)
plot(AUC~ number,data=AUC_seq_RF_16_df_09_noGC)
plot(AUC~ number,data=AUC_full_RF_16_df_09_noGC)

