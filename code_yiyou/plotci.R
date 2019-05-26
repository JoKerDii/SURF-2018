theme_set(theme_gray(base_size = 17))
m16.CI<-
matrix(nrow=5,ncol=3,NA)
colnames(m16.CI)<-c("Gen","seq","full")
m16.CI[,1]<-c(Gen_16_RF_final_06$resample$ROC)
m16.CI[,2]<-c(seq_16_RF_final_06$resample$ROC)
m16.CI[,3]<-c(full_16_RF_final_06$resample$ROC)
m16.CI<-as.data.frame(m16.CI)
S.error_m16_Gen<-sqrt(sd(m16.CI$Gen)*sd(m16.CI$Gen)/5)
S.error_m16_seq<-sqrt(sd(m16.CI$seq)*sd(m16.CI$seq)/5)
S.error_m16_full<-sqrt(sd(m16.CI$full)*sd(m16.CI$full)/5)
m16.CI_plot<-matrix(NA,nrow=3, ncol=3)
colnames(m16.CI_plot)<-c("min","mean","max")
rownames(m16.CI_plot)<-c("Gen","Seq","Full")
m16.CI_plot[1,]<-c(mean(m16.CI$Gen)-S.error_m16_Gen,mean(m16.CI$Gen),mean(m16.CI$Gen)+S.error_m16_Gen)
m16.CI_plot[2,]<-c(mean(m16.CI$seq)-S.error_m16_seq,mean(m16.CI$seq),mean(m16.CI$seq)+S.error_m16_seq)
m16.CI_plot[3,]<-c(mean(m16.CI$full)-S.error_m16_full,mean(m16.CI$full),mean(m16.CI$full)+S.error_m16_full)
m16.CI_plot<-as.data.frame(m16.CI_plot)
M16.CI_plot<-as.data.frame(m16.CI_plot$mean)
M16.CI_plot[,2]<-c("Gen","Seq","Full")
colnames(M16.CI_plot)<-c("AUC","method")
pdf("M16_06.pdf",width=7,height = 5)
ggplot(data = M16.CI_plot, aes(y=AUC,x=method))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin = m16.CI_plot$min,ymax = m16.CI_plot$max),width=0.1)+
  coord_flip()
dev.off()



theme_set(theme_gray(base_size = 17))
m16.CI<-
  matrix(nrow=5,ncol=3,NA)
colnames(m16.CI)<-c("Gen","seq","full")
m16.CI[,1]<-c(Gen_16_RF_final_06_noGC$resample$ROC)
m16.CI[,2]<-c(seq_16_RF_final_06_noGC$resample$ROC)
m16.CI[,3]<-c(full_16_RF_final_06_noGC$resample$ROC)
m16.CI<-as.data.frame(m16.CI)
S.error_m16_Gen<-sqrt(sd(m16.CI$Gen)*sd(m16.CI$Gen)/5)
S.error_m16_seq<-sqrt(sd(m16.CI$seq)*sd(m16.CI$seq)/5)
S.error_m16_full<-sqrt(sd(m16.CI$full)*sd(m16.CI$full)/5)
m16.CI_plot<-matrix(NA,nrow=3, ncol=3)
colnames(m16.CI_plot)<-c("min","mean","max")
rownames(m16.CI_plot)<-c("Gen","Seq","Full")
m16.CI_plot[1,]<-c(mean(m16.CI$Gen)-S.error_m16_Gen,mean(m16.CI$Gen),mean(m16.CI$Gen)+S.error_m16_Gen)
m16.CI_plot[2,]<-c(mean(m16.CI$seq)-S.error_m16_seq,mean(m16.CI$seq),mean(m16.CI$seq)+S.error_m16_seq)
m16.CI_plot[3,]<-c(mean(m16.CI$full)-S.error_m16_full,mean(m16.CI$full),mean(m16.CI$full)+S.error_m16_full)
m16.CI_plot<-as.data.frame(m16.CI_plot)
M16.CI_plot<-as.data.frame(m16.CI_plot$mean)
M16.CI_plot[,2]<-c("Gen","Seq","Full")
colnames(M16.CI_plot)<-c("AUC","method")
pdf("M16_06_noGC.pdf",width=7,height = 5)
ggplot(data = M16.CI_plot, aes(y=AUC,x=method))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin = m16.CI_plot$min,ymax = m16.CI_plot$max),width=0.1)+
  coord_flip()
dev.off()





writer<-matrix(NA,nrow = 3,ncol = 4)
rownames(writer)<-c("seq","Gen","full")
colnames(writer)<-c("RF","RF","RF","RF")
writer[1,]<-c(as.double(AUC_seq_RF_16_final_06),
                   as.double(AUC_seq_RF_16_final_07),
                   as.double(AUC_seq_RF_16_final_08),
                   as.double(AUC_seq_RF_16_final_09))
writer[2,]<-c(as.double(AUC_Gen_RF_16_final_06),
                   as.double(AUC_Gen_RF_16_final_07),
                   as.double(AUC_Gen_RF_16_final_08),
                   as.double(AUC_Gen_RF_16_final_09))
writer[3,]<-c(as.double(AUC_full_RF_16_final_06),
                   as.double(AUC_full_RF_16_final_07),
                   as.double(AUC_full_RF_16_final_08),
                   as.double(AUC_full_RF_16_final_09))
write.csv(writer,file = "Writer.csv")





writer_noGC<-matrix(NA,nrow = 3,ncol = 4)
rownames(writer_noGC)<-c("seq","Gen","full")
colnames(writer_noGC)<-c("RF","RF","RF","RF","RF","RF","RF","RF")
writer_noGC[1,]<-c(as.double(AUC_seq_RF_16_final_06_noGC),
                   as.double(AUC_seq_RF_16_final_07_noGC),
                   as.double(AUC_seq_RF_16_final_08_noGC),
                   as.double(AUC_seq_RF_16_final_09_noGC))
writer_noGC[2,]<-c(as.double(AUC_Gen_RF_16_final_06_noGC),
                   as.double(AUC_Gen_RF_16_final_07_noGC),
                   as.double(AUC_Gen_RF_16_final_08_noGC),
                   as.double(AUC_Gen_RF_16_final_09_noGC))
writer_noGC[3,]<-c(as.double(AUC_full_RF_16_final_06_noGC),
                   as.double(AUC_full_RF_16_final_07_noGC),
                   as.double(AUC_full_RF_16_final_08_noGC),
                   as.double(AUC_full_RF_16_final_09_noGC))
write.csv(writer_noGC,file = "Writer_noGC.csv")

M16_06_ROC<-matrix(nrow=3,ncol=3,NA)
rownames(M16_06_ROC)<-c("Gen","seq","full")
colnames(M16_06_ROC)<-c("Spec","Sens","ROC")
M16_06_ROC[1,]<-c(
  Gen_16_RF_final_06$results$Spec[which(Gen_16_RF_final_06$results$ROC==max(Gen_16_RF_final_06$results$ROC))],
  Gen_16_RF_final_06$results$Sens[which(Gen_16_RF_final_06$results$ROC==max(Gen_16_RF_final_06$results$ROC))],
  max(Gen_16_RF_final_06$results$ROC))
M16_06_ROC[2,]<-c(
  seq_16_RF_final_06$results$Spec[which(seq_16_RF_final_06$results$ROC==max(seq_16_RF_final_06$results$ROC))],
  seq_16_RF_final_06$results$Sens[which(seq_16_RF_final_06$results$ROC==max(seq_16_RF_final_06$results$ROC))],
  max(seq_16_RF_final_06$results$ROC))
M16_06_ROC[3,]<-c(
  full_16_RF_final_06$results$Spec[which(full_16_RF_final_06$results$ROC==max(full_16_RF_final_06$results$ROC))],
  full_16_RF_final_06$results$Sens[which(full_16_RF_final_06$results$ROC==max(full_16_RF_final_06$results$ROC))],
  max(full_16_RF_final_06$results$ROC))
write.csv(M16_06_ROC,"M16_06_ROC.csv")




M16_06_noGC_ROC<-matrix(nrow=3,ncol=3,NA)
rownames(M16_06_noGC_ROC)<-c("Gen","seq","full")
colnames(M16_06_noGC_ROC)<-c("Spec","Sens","ROC")
M16_06_noGC_ROC[1,]<-c(
  Gen_16_RF_final_06_noGC$results$Spec[which(Gen_16_RF_final_06_noGC$results$ROC==max(Gen_16_RF_final_06_noGC$results$ROC))],
  Gen_16_RF_final_06_noGC$results$Sens[which(Gen_16_RF_final_06_noGC$results$ROC==max(Gen_16_RF_final_06_noGC$results$ROC))],
  max(Gen_16_RF_final_06_noGC$results$ROC))
M16_06_noGC_ROC[2,]<-c(
  seq_16_RF_final_06_noGC$results$Spec[which(seq_16_RF_final_06_noGC$results$ROC==max(seq_16_RF_final_06_noGC$results$ROC))],
  seq_16_RF_final_06_noGC$results$Sens[which(seq_16_RF_final_06_noGC$results$ROC==max(seq_16_RF_final_06_noGC$results$ROC))],
  max(seq_16_RF_final_06_noGC$results$ROC))
M16_06_noGC_ROC[3,]<-c(
  full_16_RF_final_06_noGC$results$Spec[which(full_16_RF_final_06_noGC$results$ROC==max(full_16_RF_final_06_noGC$results$ROC))],
  full_16_RF_final_06_noGC$results$Sens[which(full_16_RF_final_06_noGC$results$ROC==max(full_16_RF_final_06_noGC$results$ROC))],
  max(full_16_RF_final_06_noGC$results$ROC))
write.csv(M16_06_noGC_ROC,"M16_06_noGC_ROC.csv")

