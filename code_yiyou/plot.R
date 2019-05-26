theme_set(theme_gray(base_size = 17))

  ITcRF3<-matrix(NA,nrow=12, ncol=5)
colnames(ITcRF3)<-c("AUC","method","site_probability","min","max")
ITcRF3[,1]<-c(as.double(ROC_Gen_RF_3_final_06@y.values),as.double(ROC_Gen_RF_3_final_07@y.values),  
as.double(ROC_Gen_RF_3_final_08@y.values),as.double(ROC_Gen_RF_3_final_09@y.values),
as.double(ROC_seq_RF_3_final_06@y.values),as.double(ROC_seq_RF_3_final_07@y.values),
as.double(ROC_seq_RF_3_final_08@y.values),as.double(ROC_seq_RF_3_final_09@y.values),
as.double(ROC_full_RF_3_final_06@y.values),as.double(ROC_full_RF_3_final_07@y.values),
as.double(ROC_full_RF_3_final_08@y.values),as.double(ROC_full_RF_3_final_09@y.values))
ITcRF3[,2]<-c("bio","bio","bio","bio","seq","seq","seq","seq","full","full","full","full")
ITcRF3[,3]<-c("0.6","0.7","0.8","0.9","0.6","0.7","0.8","0.9","0.6","0.7","0.8","0.9")
ITcRF3[,4]<-c(c(as.double(ROC_Gen_RF_3_final_06@y.values)-sqrt(var(roc(Gen_testing3_06$M3,Gen_testing3_06$RFpredfin))),
                 as.double(ROC_Gen_RF_3_final_07@y.values)-sqrt(var(roc(Gen_testing3_07$M3,Gen_testing3_07$RFpredfin))),  
                 as.double(ROC_Gen_RF_3_final_08@y.values)-sqrt(var(roc(Gen_testing3_08$M3,Gen_testing3_08$RFpredfin))),
                 as.double(ROC_Gen_RF_3_final_09@y.values)-sqrt(var(roc(Gen_testing3_09$M3,Gen_testing3_09$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_06@y.values)-sqrt(var(roc(seq_testing3_06$M3,seq_testing3_06$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_07@y.values)-sqrt(var(roc(seq_testing3_07$M3,seq_testing3_07$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_08@y.values)-sqrt(var(roc(seq_testing3_08$M3,seq_testing3_08$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_09@y.values)-sqrt(var(roc(seq_testing3_09$M3,seq_testing3_09$RFpredfin))),
                 as.double(ROC_full_RF_3_final_06@y.values)-sqrt(var(roc(full_testing3_06$M3,full_testing3_06$RFpredfin))),
                 as.double(ROC_full_RF_3_final_07@y.values)-sqrt(var(roc(full_testing3_07$M3,full_testing3_07$RFpredfin))),
                 as.double(ROC_full_RF_3_final_08@y.values)-sqrt(var(roc(full_testing3_08$M3,full_testing3_08$RFpredfin))),
                 as.double(ROC_full_RF_3_final_09@y.values)-sqrt(var(roc(full_testing3_09$M3,full_testing3_09$RFpredfin)))))
ITcRF3[,5]<-c(c(as.double(ROC_Gen_RF_3_final_06@y.values)+sqrt(var(roc(Gen_testing3_06$M3,Gen_testing3_06$RFpredfin))),
                 as.double(ROC_Gen_RF_3_final_07@y.values)+sqrt(var(roc(Gen_testing3_07$M3,Gen_testing3_07$RFpredfin))),  
                 as.double(ROC_Gen_RF_3_final_08@y.values)+sqrt(var(roc(Gen_testing3_08$M3,Gen_testing3_08$RFpredfin))),
                 as.double(ROC_Gen_RF_3_final_09@y.values)+sqrt(var(roc(Gen_testing3_09$M3,Gen_testing3_09$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_06@y.values)+sqrt(var(roc(seq_testing3_06$M3,seq_testing3_06$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_07@y.values)+sqrt(var(roc(seq_testing3_07$M3,seq_testing3_07$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_08@y.values)+sqrt(var(roc(seq_testing3_08$M3,seq_testing3_08$RFpredfin))),
                 as.double(ROC_seq_RF_3_final_09@y.values)+sqrt(var(roc(seq_testing3_09$M3,seq_testing3_09$RFpredfin))),
                 as.double(ROC_full_RF_3_final_06@y.values)+sqrt(var(roc(full_testing3_06$M3,full_testing3_06$RFpredfin))),
                 as.double(ROC_full_RF_3_final_07@y.values)+sqrt(var(roc(full_testing3_07$M3,full_testing3_07$RFpredfin))),
                 as.double(ROC_full_RF_3_final_08@y.values)+sqrt(var(roc(full_testing3_08$M3,full_testing3_08$RFpredfin))),
                 as.double(ROC_full_RF_3_final_09@y.values)+sqrt(var(roc(full_testing3_09$M3,full_testing3_09$RFpredfin)))))
ITcRF3<-as.data.frame(ITcRF3)
ITcRF3$AUC<-unfactor(ITcRF3$AUC)
ITcRF3$min<-unfactor(ITcRF3$min)
ITcRF3$max<-unfactor(ITcRF3$max)
pdf("M3.pdf",width=7,height = 5) 
ggplot(data = ITcRF3, aes(x=site_probability,y=AUC,color = method))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin = ITcRF3$min,ymax = ITcRF3$max),width=0.1)+
  scale_y_continuous(name="AUC")+
labs(x="site probability",title = "AUC score for different models (METTL3)")
dev.off()


ITcRF14<-matrix(NA,nrow=12, ncol=5)
colnames(ITcRF14)<-c("AUC","method","site_probability","min","max")
ITcRF14[,1]<-c(as.double(ROC_Gen_RF_14_final_06@y.values),as.double(ROC_Gen_RF_14_final_07@y.values),  
               as.double(ROC_Gen_RF_14_final_08@y.values),as.double(ROC_Gen_RF_14_final_09@y.values),
               as.double(ROC_seq_RF_14_final_06@y.values),as.double(ROC_seq_RF_14_final_07@y.values),
               as.double(ROC_seq_RF_14_final_08@y.values),as.double(ROC_seq_RF_14_final_09@y.values),
               as.double(ROC_full_RF_14_final_06@y.values),as.double(ROC_full_RF_14_final_07@y.values),
               as.double(ROC_full_RF_14_final_08@y.values),as.double(ROC_full_RF_14_final_09@y.values))
ITcRF14[,2]<-c("bio","bio","bio","bio","seq","seq","seq","seq","full","full","full","full")
ITcRF14[,3]<-c("0.6","0.7","0.8","0.9","0.6","0.7","0.8","0.9","0.6","0.7","0.8","0.9")
ITcRF14[,4]<-c(c(as.double(ROC_Gen_RF_14_final_06@y.values)-sqrt(var(roc(Gen_testing14_06$M14,Gen_testing14_06$RFpredfin))),
                 as.double(ROC_Gen_RF_14_final_07@y.values)-sqrt(var(roc(Gen_testing14_07$M14,Gen_testing14_07$RFpredfin))),  
                 as.double(ROC_Gen_RF_14_final_08@y.values)-sqrt(var(roc(Gen_testing14_08$M14,Gen_testing14_08$RFpredfin))),
                 as.double(ROC_Gen_RF_14_final_09@y.values)-sqrt(var(roc(Gen_testing14_09$M14,Gen_testing14_09$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_06@y.values)-sqrt(var(roc(seq_testing14_06$M14,seq_testing14_06$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_07@y.values)-sqrt(var(roc(seq_testing14_07$M14,seq_testing14_07$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_08@y.values)-sqrt(var(roc(seq_testing14_08$M14,seq_testing14_08$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_09@y.values)-sqrt(var(roc(seq_testing14_09$M14,seq_testing14_09$RFpredfin))),
                 as.double(ROC_full_RF_14_final_06@y.values)-sqrt(var(roc(full_testing14_06$M14,full_testing14_06$RFpredfin))),
                 as.double(ROC_full_RF_14_final_07@y.values)-sqrt(var(roc(full_testing14_07$M14,full_testing14_07$RFpredfin))),
                 as.double(ROC_full_RF_14_final_08@y.values)-sqrt(var(roc(full_testing14_08$M14,full_testing14_08$RFpredfin))),
                 as.double(ROC_full_RF_14_final_09@y.values)-sqrt(var(roc(full_testing14_09$M14,full_testing14_09$RFpredfin)))))
ITcRF14[,5]<-c(c(as.double(ROC_Gen_RF_14_final_06@y.values)+sqrt(var(roc(Gen_testing14_06$M14,Gen_testing14_06$RFpredfin))),
                 as.double(ROC_Gen_RF_14_final_07@y.values)+sqrt(var(roc(Gen_testing14_07$M14,Gen_testing14_07$RFpredfin))),  
                 as.double(ROC_Gen_RF_14_final_08@y.values)+sqrt(var(roc(Gen_testing14_08$M14,Gen_testing14_08$RFpredfin))),
                 as.double(ROC_Gen_RF_14_final_09@y.values)+sqrt(var(roc(Gen_testing14_09$M14,Gen_testing14_09$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_06@y.values)+sqrt(var(roc(seq_testing14_06$M14,seq_testing14_06$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_07@y.values)+sqrt(var(roc(seq_testing14_07$M14,seq_testing14_07$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_08@y.values)+sqrt(var(roc(seq_testing14_08$M14,seq_testing14_08$RFpredfin))),
                 as.double(ROC_seq_RF_14_final_09@y.values)+sqrt(var(roc(seq_testing14_09$M14,seq_testing14_09$RFpredfin))),
                 as.double(ROC_full_RF_14_final_06@y.values)+sqrt(var(roc(full_testing14_06$M14,full_testing14_06$RFpredfin))),
                 as.double(ROC_full_RF_14_final_07@y.values)+sqrt(var(roc(full_testing14_07$M14,full_testing14_07$RFpredfin))),
                 as.double(ROC_full_RF_14_final_08@y.values)+sqrt(var(roc(full_testing14_08$M14,full_testing14_08$RFpredfin))),
                 as.double(ROC_full_RF_14_final_09@y.values)+sqrt(var(roc(full_testing14_09$M14,full_testing14_09$RFpredfin)))))
ITcRF14<-as.data.frame(ITcRF14)
ITcRF14$AUC<-unfactor(ITcRF14$AUC)
ITcRF14$min<-unfactor(ITcRF14$min)
ITcRF14$max<-unfactor(ITcRF14$max)
pdf("M14.pdf",width=7,height = 5) 
ggplot(data = ITcRF14, aes(x=site_probability,y=AUC,color = method))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin = ITcRF14$min,ymax = ITcRF14$max),width=0.1)+
  scale_y_continuous(name="AUC")+
  labs(x="site probability",title = "AUC score for different models (METTL14)")
dev.off()







ITcRF16<-matrix(NA,nrow=12, ncol=5)
colnames(ITcRF16)<-c("AUC","method","site_probability","min","max")
ITcRF16[,1]<-c(as.double(ROC_Gen_RF_16_final_06@y.values),as.double(ROC_Gen_RF_16_final_07@y.values),  
               as.double(ROC_Gen_RF_16_final_08@y.values),as.double(ROC_Gen_RF_16_final_09@y.values),
               as.double(ROC_seq_RF_16_final_06@y.values),as.double(ROC_seq_RF_16_final_07@y.values),
               as.double(ROC_seq_RF_16_final_08@y.values),as.double(ROC_seq_RF_16_final_09@y.values),
               as.double(ROC_full_RF_16_final_06@y.values),as.double(ROC_full_RF_16_final_07@y.values),
               as.double(ROC_full_RF_16_final_08@y.values),as.double(ROC_full_RF_16_final_09@y.values))
ITcRF16[,2]<-c("bio","bio","bio","bio","seq","seq","seq","seq","full","full","full","full")
ITcRF16[,3]<-c("0.6","0.7","0.8","0.9","0.6","0.7","0.8","0.9","0.6","0.7","0.8","0.9")
ITcRF16[,4]<-c(c(as.double(ROC_Gen_RF_16_final_06@y.values)-sqrt(var(roc(Gen_testing16_06$M16,Gen_testing16_06$RFpredfin))),
                 as.double(ROC_Gen_RF_16_final_07@y.values)-sqrt(var(roc(Gen_testing16_07$M16,Gen_testing16_07$RFpredfin))),  
                 as.double(ROC_Gen_RF_16_final_08@y.values)-sqrt(var(roc(Gen_testing16_08$M16,Gen_testing16_08$RFpredfin))),
                 as.double(ROC_Gen_RF_16_final_09@y.values)-sqrt(var(roc(Gen_testing16_09$M16,Gen_testing16_09$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_06@y.values)-sqrt(var(roc(seq_testing16_06$M16,seq_testing16_06$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_07@y.values)-sqrt(var(roc(seq_testing16_07$M16,seq_testing16_07$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_08@y.values)-sqrt(var(roc(seq_testing16_08$M16,seq_testing16_08$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_09@y.values)-sqrt(var(roc(seq_testing16_09$M16,seq_testing16_09$RFpredfin))),
                 as.double(ROC_full_RF_16_final_06@y.values)-sqrt(var(roc(full_testing16_06$M16,full_testing16_06$RFpredfin))),
                 as.double(ROC_full_RF_16_final_07@y.values)-sqrt(var(roc(full_testing16_07$M16,full_testing16_07$RFpredfin))),
                 as.double(ROC_full_RF_16_final_08@y.values)-sqrt(var(roc(full_testing16_08$M16,full_testing16_08$RFpredfin))),
                 as.double(ROC_full_RF_16_final_09@y.values)-sqrt(var(roc(full_testing16_09$M16,full_testing16_09$RFpredfin)))))
ITcRF16[,5]<-c(c(as.double(ROC_Gen_RF_16_final_06@y.values)+sqrt(var(roc(Gen_testing16_06$M16,Gen_testing16_06$RFpredfin))),
                 as.double(ROC_Gen_RF_16_final_07@y.values)+sqrt(var(roc(Gen_testing16_07$M16,Gen_testing16_07$RFpredfin))),  
                 as.double(ROC_Gen_RF_16_final_08@y.values)+sqrt(var(roc(Gen_testing16_08$M16,Gen_testing16_08$RFpredfin))),
                 as.double(ROC_Gen_RF_16_final_09@y.values)+sqrt(var(roc(Gen_testing16_09$M16,Gen_testing16_09$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_06@y.values)+sqrt(var(roc(seq_testing16_06$M16,seq_testing16_06$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_07@y.values)+sqrt(var(roc(seq_testing16_07$M16,seq_testing16_07$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_08@y.values)+sqrt(var(roc(seq_testing16_08$M16,seq_testing16_08$RFpredfin))),
                 as.double(ROC_seq_RF_16_final_09@y.values)+sqrt(var(roc(seq_testing16_09$M16,seq_testing16_09$RFpredfin))),
                 as.double(ROC_full_RF_16_final_06@y.values)+sqrt(var(roc(full_testing16_06$M16,full_testing16_06$RFpredfin))),
                 as.double(ROC_full_RF_16_final_07@y.values)+sqrt(var(roc(full_testing16_07$M16,full_testing16_07$RFpredfin))),
                 as.double(ROC_full_RF_16_final_08@y.values)+sqrt(var(roc(full_testing16_08$M16,full_testing16_08$RFpredfin))),
                 as.double(ROC_full_RF_16_final_09@y.values)+sqrt(var(roc(full_testing16_09$M16,full_testing16_09$RFpredfin)))))
ITcRF16<-as.data.frame(ITcRF16)
ITcRF16$AUC<-unfactor(ITcRF16$AUC)
ITcRF16$min<-unfactor(ITcRF16$min)
ITcRF16$max<-unfactor(ITcRF16$max)
pdf("M16.pdf",width=7,height = 5)
ggplot(data = ITcRF16, aes(x=site_probability,y=AUC,color = method))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin = ITcRF16$min,ymax = ITcRF16$max),width=0.1)+
  scale_y_continuous(name="AUC")+
  labs(x="site probability",title = "AUC score for different models (METTL16)")
dev.off()





AUC_Gen_RF_3_df_09<-data.frame(AUC=AUC_Gen_RF_3v_09,number=c(1:62))
AUC_seq_RF_3_df_09<-data.frame(AUC=AUC_seq_RF_3v_09,number=c(1:155))
AUC_full_RF_3_df_09<-data.frame(AUC=AUC_full_RF_3v_09,number=c(1:217))
AUC_Gen_RF_14_df_09<-data.frame(AUC=AUC_Gen_RF_14v_09,number=c(1:62))
AUC_seq_RF_14_df_09<-data.frame(AUC=AUC_seq_RF_14v_09,number=c(1:155))
AUC_full_RF_14_df_09<-data.frame(AUC=AUC_full_RF_14v_09,number=c(1:217))
AUC_Gen_RF_16_df_09<-data.frame(AUC=AUC_Gen_RF_16v_09,number=c(1:62))
AUC_seq_RF_16_df_09<-data.frame(AUC=AUC_seq_RF_16v_09,number=c(1:155))
AUC_full_RF_16_df_09<-data.frame(AUC=AUC_full_RF_16v_09,number=c(1:217))

AUC_Gen_RF_3_df_08<-data.frame(AUC=AUC_Gen_RF_3v_08,number=c(1:62))
AUC_seq_RF_3_df_08<-data.frame(AUC=AUC_seq_RF_3v_08,number=c(1:155))
AUC_full_RF_3_df_08<-data.frame(AUC=AUC_full_RF_3v_08,number=c(1:217))
AUC_Gen_RF_14_df_08<-data.frame(AUC=AUC_Gen_RF_14v_08,number=c(1:62))
AUC_seq_RF_14_df_08<-data.frame(AUC=AUC_seq_RF_14v_08,number=c(1:155))
AUC_full_RF_14_df_08<-data.frame(AUC=AUC_full_RF_14v_08,number=c(1:217))
AUC_Gen_RF_16_df_08<-data.frame(AUC=AUC_Gen_RF_16v_08,number=c(1:62))
AUC_seq_RF_16_df_08<-data.frame(AUC=AUC_seq_RF_16v_08,number=c(1:155))
AUC_full_RF_16_df_08<-data.frame(AUC=AUC_full_RF_16v_08,number=c(1:217))


AUC_Gen_RF_3_df_07<-data.frame(AUC=AUC_Gen_RF_3v_07,number=c(1:62))
AUC_seq_RF_3_df_07<-data.frame(AUC=AUC_seq_RF_3v_07,number=c(1:155))
AUC_full_RF_3_df_07<-data.frame(AUC=AUC_full_RF_3v_07,number=c(1:217))
AUC_Gen_RF_14_df_07<-data.frame(AUC=AUC_Gen_RF_14v_07,number=c(1:62))
AUC_seq_RF_14_df_07<-data.frame(AUC=AUC_seq_RF_14v_07,number=c(1:155))
AUC_full_RF_14_df_07<-data.frame(AUC=AUC_full_RF_14v_07,number=c(1:217))
AUC_Gen_RF_16_df_07<-data.frame(AUC=AUC_Gen_RF_16v_07,number=c(1:62))
AUC_seq_RF_16_df_07<-data.frame(AUC=AUC_seq_RF_16v_07,number=c(1:155))
AUC_full_RF_16_df_07<-data.frame(AUC=AUC_full_RF_16v_07,number=c(1:217))

AUC_Gen_RF_3_df_06<-data.frame(AUC=AUC_Gen_RF_3v_06,number=c(1:62))
AUC_seq_RF_3_df_06<-data.frame(AUC=AUC_seq_RF_3v_06,number=c(1:155))
AUC_full_RF_3_df_06<-data.frame(AUC=AUC_full_RF_3v_06,number=c(1:217))
AUC_Gen_RF_14_df_06<-data.frame(AUC=AUC_Gen_RF_14v_06,number=c(1:62))
AUC_seq_RF_14_df_06<-data.frame(AUC=AUC_seq_RF_14v_06,number=c(1:155))
AUC_full_RF_14_df_06<-data.frame(AUC=AUC_full_RF_14v_06,number=c(1:217))
AUC_Gen_RF_16_df_06<-data.frame(AUC=AUC_Gen_RF_16v_06,number=c(1:62))
AUC_seq_RF_16_df_06<-data.frame(AUC=AUC_seq_RF_16v_06,number=c(1:155))
AUC_full_RF_16_df_06<-data.frame(AUC=AUC_full_RF_16v_06,number=c(1:217))



pdf("AUC_scat_full_06_3.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06)))
a[which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))]<-a[which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))]-8
a[1:which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))-1]<-a[1:which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))-1]-7
a=a+8
ggplot(data=AUC_full_RF_3_df_06, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL3, cut-off>0.6)")+
  theme(plot.title = element_text(size=15))
dev.off()



pdf("AUC_scat_full_06_14.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06)))
a[which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))]<-a[which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))]-8
a[1:which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))-1]<-a[1:which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))-1]-7
a=a+8
ggplot(data=AUC_full_RF_14_df_06, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL14, cut-off>0.6)")+
  theme(plot.title = element_text(size=15))
dev.off()





pdf("AUC_scat_full_06_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06)))
a[which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))]<-a[which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))]-8
a[1:which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))-1]<-a[1:which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))-1]-7
a=a+8
ggplot(data=AUC_full_RF_16_df_06, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.6)")+
  theme(plot.title = element_text(size=15))
dev.off()


pdf("importance score full_06_3.pdf",width=7,height = 5) 

full_imp3_06$variable <-rownames(full_imp3_06)
full_imp3_06$variable <- factor(full_imp3_06$variable,levels = full_imp3_06$variable)

a<-as.integer(as.logical(AUC_full_RF_3v_06[1:15]==max(AUC_full_RF_3v_06)))
a[1:which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))-1]<-a[1:which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))-1]+1
a[(1+which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))):15]<-a[(1+which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))):15]+8

ggplot(data=data.frame(Overall=full_imp3_06$Overall[1:15],variable=full_imp3_06$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL3, cut-off>0.6)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()




pdf("importance score full_06_14.pdf",width=7,height = 5) 

full_imp14_06$variable <-rownames(full_imp14_06)
full_imp14_06$variable <- factor(full_imp14_06$variable,levels = full_imp14_06$variable)

a<-as.integer(as.logical(AUC_full_RF_14v_06[1:15]==max(AUC_full_RF_14v_06)))
a[1:which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))-1]<-a[1:which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))-1]+1
a[(1+which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))):15]<-a[(1+which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))):15]+8

ggplot(data=data.frame(Overall=full_imp14_06$Overall[1:15],variable=full_imp14_06$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL14, cut-off>0.6)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("importance score full_06_16.pdf",width=7,height = 5) 

full_imp16_06$variable <-rownames(full_imp16_06)
full_imp16_06$variable <- factor(full_imp16_06$variable,levels = full_imp16_06$variable)

a<-as.integer(as.logical(AUC_full_RF_16v_06[1:15]==max(AUC_full_RF_16v_06)))
a[1:which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))-1]<-a[1:which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))-1]+1
a[(1+which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))):15]<-a[(1+which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))):15]+8

ggplot(data=data.frame(Overall=full_imp16_06$Overall[1:15],variable=full_imp16_06$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16, cut-off>0.6)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()

























methyltransferase<- read.csv("methyltransferase.csv",stringsAsFactors = F)[,c("Term","PValue","methyltransferase")]
for (i in 1:nrow(methyltransferase) ){
  index <- gregexpr("~",methyltransferase[i,])[[1]][1]
  methyltransferase[i,1] <- substr(methyltransferase[i,1],(index+1),nchar(c(methyltransferase[i,1])))
}

methyltransferase<-methyltransferase[order(methyltransferase$methyltransferase),,drop=FALSE]
methyltransferase$Term <- factor(methyltransferase$Term,levels = methyltransferase$Term)
pdf("GO_plot.pdf",width=7,height = 5) 
ggplot(data = methyltransferase, aes(x=methyltransferase,y=Term,color = PValue))+
  geom_point(aes(size=-log(PValue)))+
  labs(title="Functions of METLL3/14/16")+
  scale_size_continuous(range = c(2,7))+
  theme_gray(base_size = 11)
dev.off()



# pdf("AUC plot_06.pdf",width = 12,height = 3)
# ggplot(data=csv09[13:18,], aes(x=rownames(csv06[13:18,]), y=AUC))+
#   geom_bar(stat = "identity",fill=c("white","gray","white","grey","white","grey"))+
#   geom_text(aes(label=csv09$AUC[13:18]))+
#   labs(x="model",title = "AUC score for sites with probablity>0.9")+
#   theme(axis.text=element_text(size=8))+
#   coord_flip()
# dev.off()
# 

summary(resamples(list(RF16= full_16_RF_final_06,RF14= full_14_RF_final_06,RF3= full_3_RF_final_06)))

sqrt(var(roc(full_testing14_09$M14,full_testing14_09$RFpredfin)))


which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))
which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))
which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))

which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))
which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))
which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))

which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))
which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))
which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))

which(AUC_full_RF_16v_06==max(AUC_full_RF_16v_06))
which(AUC_full_RF_14v_06==max(AUC_full_RF_14v_06))
which(AUC_full_RF_3v_06==max(AUC_full_RF_3v_06))







pdf("AUC_scat_full_07_3.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07)))
a[which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))]<-a[which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))]-8
a[1:which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))-1]<-a[1:which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))-1]-7
a=a+8
ggplot(data=AUC_full_RF_3_df_07, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL3, cut-off>0.7)")+
  theme(plot.title = element_text(size=15))
dev.off()



pdf("AUC_scat_full_07_14.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07)))
a[which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))]<-a[which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))]-8
a[1:which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))-1]<-a[1:which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))-1]-7
a=a+8
ggplot(data=AUC_full_RF_14_df_07, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL14, cut-off>0.7)")+
  theme(plot.title = element_text(size=15))
dev.off()





pdf("AUC_scat_full_07_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07)))
a[which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))]<-a[which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))]-8
a[1:which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))-1]<-a[1:which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))-1]-7
a=a+8
ggplot(data=AUC_full_RF_16_df_07, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.7)")+
  theme(plot.title = element_text(size=15))
dev.off()


pdf("importance score full_07_3.pdf",width=7,height = 5) 

full_imp3_07$variable <-rownames(full_imp3_07)
full_imp3_07$variable <- factor(full_imp3_07$variable,levels = full_imp3_07$variable)

a<-as.integer(as.logical(AUC_full_RF_3v_07[1:15]==max(AUC_full_RF_3v_07)))
a[1:which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))-1]<-a[1:which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))-1]+1
a[(1+which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))):15]<-a[(1+which(AUC_full_RF_3v_07==max(AUC_full_RF_3v_07))):15]+8

ggplot(data=data.frame(Overall=full_imp3_07$Overall[1:15],variable=full_imp3_07$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL3, cut-off>0.7)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()




pdf("importance score full_07_14.pdf",width=7,height = 5) 

full_imp14_07$variable <-rownames(full_imp14_07)
full_imp14_07$variable <- factor(full_imp14_07$variable,levels = full_imp14_07$variable)

a<-as.integer(as.logical(AUC_full_RF_14v_07[1:15]==max(AUC_full_RF_14v_07)))
a[1:which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))-1]<-a[1:which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))-1]+1
a[(1+which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))):15]<-a[(1+which(AUC_full_RF_14v_07==max(AUC_full_RF_14v_07))):15]+8

ggplot(data=data.frame(Overall=full_imp14_07$Overall[1:15],variable=full_imp14_07$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL14, cut-off>0.7)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("importance score full_07_16.pdf",width=7,height = 5) 

full_imp16_07$variable <-rownames(full_imp16_07)
full_imp16_07$variable <- factor(full_imp16_07$variable,levels = full_imp16_07$variable)

a<-as.integer(as.logical(AUC_full_RF_16v_07[1:15]==max(AUC_full_RF_16v_07)))
a[1:which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))-1]<-a[1:which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))-1]+1
a[(1+which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))):15]<-a[(1+which(AUC_full_RF_16v_07==max(AUC_full_RF_16v_07))):15]+8

ggplot(data=data.frame(Overall=full_imp16_07$Overall[1:15],variable=full_imp16_07$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16, cut-off>0.7)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()


























pdf("AUC_scat_full_08_3.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08)))
a[which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))]<-a[which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))]-8
a[1:which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))-1]<-a[1:which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))-1]-7
a=a+8
ggplot(data=AUC_full_RF_3_df_08, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL3, cut-off>0.8)")+
  theme(plot.title = element_text(size=15))
dev.off()



pdf("AUC_scat_full_08_14.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08)))
a[which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))]<-a[which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))]-8
a[1:which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))-1]<-a[1:which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))-1]-7
a=a+8
ggplot(data=AUC_full_RF_14_df_08, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL14, cut-off>0.8)")+
  theme(plot.title = element_text(size=15))
dev.off()





pdf("AUC_scat_full_08_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08)))
a[which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))]<-a[which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))]-8
a[1:which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))-1]<-a[1:which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))-1]-7
a=a+8
ggplot(data=AUC_full_RF_16_df_08, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.8)")+
  theme(plot.title = element_text(size=15))
dev.off()


pdf("importance score full_08_3.pdf",width=7,height = 5) 

full_imp3_08$variable <-rownames(full_imp3_08)
full_imp3_08$variable <- factor(full_imp3_08$variable,levels = full_imp3_08$variable)

a<-as.integer(as.logical(AUC_full_RF_3v_08[1:15]==max(AUC_full_RF_3v_08)))
a[1:which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))-1]<-a[1:which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))-1]+1
a[(1+which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))):15]<-a[(1+which(AUC_full_RF_3v_08==max(AUC_full_RF_3v_08))):15]+8

ggplot(data=data.frame(Overall=full_imp3_08$Overall[1:15],variable=full_imp3_08$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL3, cut-off>0.8)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()




pdf("importance score full_08_14.pdf",width=7,height = 5) 

full_imp14_08$variable <-rownames(full_imp14_08)
full_imp14_08$variable <- factor(full_imp14_08$variable,levels = full_imp14_08$variable)

a<-as.integer(as.logical(AUC_full_RF_14v_08[1:15]==max(AUC_full_RF_14v_08)))
a[1:which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))-1]<-a[1:which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))-1]+1
a[(1+which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))):15]<-a[(1+which(AUC_full_RF_14v_08==max(AUC_full_RF_14v_08))):15]+8

ggplot(data=data.frame(Overall=full_imp14_08$Overall[1:15],variable=full_imp14_08$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL14, cut-off>0.8)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("importance score full_08_16.pdf",width=7,height = 5) 

full_imp16_08$variable <-rownames(full_imp16_08)
full_imp16_08$variable <- factor(full_imp16_08$variable,levels = full_imp16_08$variable)

a<-as.integer(as.logical(AUC_full_RF_16v_08[1:15]==max(AUC_full_RF_16v_08)))
a[1:which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))-1]<-a[1:which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))-1]+1
a[(1+which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))):15]<-a[(1+which(AUC_full_RF_16v_08==max(AUC_full_RF_16v_08))):15]+8

ggplot(data=data.frame(Overall=full_imp16_08$Overall[1:15],variable=full_imp16_08$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16, cut-off>0.8)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()




























pdf("AUC_scat_full_09_3.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09)))
a[which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))]<-a[which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))]-8
a[1:which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))-1]<-a[1:which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))-1]-7
a=a+8
ggplot(data=AUC_full_RF_3_df_09, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL3, cut-off>0.9)")+
  theme(plot.title = element_text(size=15))
dev.off()



pdf("AUC_scat_full_09_14.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09)))
a[which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))]<-a[which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))]-8
a[1:which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))-1]<-a[1:which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))-1]-7
a=a+8
ggplot(data=AUC_full_RF_14_df_09, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL14, cut-off>0.9)")+
  theme(plot.title = element_text(size=15))
dev.off()





pdf("AUC_scat_full_09_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09)))
a[which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))]<-a[which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))]-8
a[1:which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))-1]<-a[1:which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))-1]-7
a=a+8
ggplot(data=AUC_full_RF_16_df_09, aes(x=log(number), y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.9)")+
  theme(plot.title = element_text(size=15))
dev.off()


pdf("importance score full_09_3.pdf",width=7,height = 5) 

full_imp3_09$variable <-rownames(full_imp3_09)
full_imp3_09$variable <- factor(full_imp3_09$variable,levels = full_imp3_09$variable)

a<-as.integer(as.logical(AUC_full_RF_3v_09[1:15]==max(AUC_full_RF_3v_09)))
a[1:which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))-1]<-a[1:which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))-1]+1
a[(1+which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))):15]<-a[(1+which(AUC_full_RF_3v_09==max(AUC_full_RF_3v_09))):15]+8

ggplot(data=data.frame(Overall=full_imp3_09$Overall[1:15],variable=full_imp3_09$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL3, cut-off>0.9)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()




pdf("importance score full_09_14.pdf",width=7,height = 5) 

full_imp14_09$variable <-rownames(full_imp14_09)
full_imp14_09$variable <- factor(full_imp14_09$variable,levels = full_imp14_09$variable)

a<-as.integer(as.logical(AUC_full_RF_14v_09[1:15]==max(AUC_full_RF_14v_09)))
a[1:which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))-1]<-a[1:which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))-1]+1
a[(1+which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))):15]<-a[(1+which(AUC_full_RF_14v_09==max(AUC_full_RF_14v_09))):15]+8

ggplot(data=data.frame(Overall=full_imp14_09$Overall[1:15],variable=full_imp14_09$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL14, cut-off>0.9)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("importance score full_09_16.pdf",width=7,height = 5) 

full_imp16_09$variable <-rownames(full_imp16_09)
full_imp16_09$variable <- factor(full_imp16_09$variable,levels = full_imp16_09$variable)

a<-as.integer(as.logical(AUC_full_RF_16v_09[1:15]==max(AUC_full_RF_16v_09)))
a[1:which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))-1]<-a[1:which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))-1]+1
a[(1+which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))):15]<-a[(1+which(AUC_full_RF_16v_09==max(AUC_full_RF_16v_09))):15]+8

ggplot(data=data.frame(Overall=full_imp16_09$Overall[1:15],variable=full_imp16_09$variable[1:15]), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16,cut-off>0.9)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()












resamps <- resamples(list(Gen=Gen_16_RF_final_06,seq=seq_16_RF_final_06,full=full_16_RF_final_06))
summary(resamps)
resamps <- resamples(list(Gen=Gen_16_RF_final_06_noGC,seq=seq_16_RF_final_06_noGC,full=full_16_RF_final_06_noGC))
summary(resamps)









pdf("importance score Gen_06_16.pdf",width=7,height = 5) 

Gen_imp16_06$variable <-rownames(Gen_imp16_06)
Gen_imp16_06$variable <- factor(Gen_imp16_06$variable,levels = Gen_imp16_06$variable)

a<-as.integer(as.logical(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06)))
a[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))-1]<-a[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))-1]+1
a[(1+which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))):50]<-a[(1+which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))):50]+8

ggplot(data=data.frame(Overall=Gen_imp16_06$Overall,variable=Gen_imp16_06$variable), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16, cut-off>0.6)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("AUC_scat_Gen_06_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06)))
a[which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))]<-a[which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))]-8
a[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))-1]<-a[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))-1]-7
a=a+8
ggplot(data=AUC_Gen_RF_16_df_06, aes(x=number, y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.6)")+
  theme(plot.title = element_text(size=15))
dev.off()







pdf("importance score Gen_06_noGC_16.pdf",width=7,height = 5) 

Gen_imp16_06_noGC$variable <-rownames(Gen_imp16_06_noGC)
Gen_imp16_06_noGC$variable <- factor(Gen_imp16_06_noGC$variable,levels = Gen_imp16_06_noGC$variable)

a<-as.integer(as.logical(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC)))
a[1:which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))-1]<-a[1:which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))-1]+1
a[(1+which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))):47]<-a[(1+which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))):47]+1

ggplot(data=data.frame(Overall=Gen_imp16_06_noGC$Overall,variable=Gen_imp16_06_noGC$variable), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16, cut-off>0.6)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("AUC_scat_Gen_06_noGC_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC)))
a[which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))]<-a[which(AUC_Gen_RF_16v_06_noGC==max(AUC_Gen_RF_16v_06_noGC))]-1
a<-a-7
a=a+8
ggplot(data=AUC_Gen_RF_16_df_06_noGC, aes(x=number, y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.6)")+
  theme(plot.title = element_text(size=15))
dev.off()





pdf("importance score Gen_06_16.pdf",width=7,height = 5) 

Gen_imp16_06$variable <-rownames(Gen_imp16_06)
Gen_imp16_06$variable <- factor(Gen_imp16_06$variable,levels = Gen_imp16_06$variable)
a<-0
a<-as.integer(as.logical(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06)))
a[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))-1]<-a[1:which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))-1]+1
a[(1+which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))):50]<-a[(1+which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))):50]+1

ggplot(data=data.frame(Overall=Gen_imp16_06$Overall,variable=Gen_imp16_06$variable), aes(x=variable, y=Overall))+
  geom_bar(stat = "identity",fill=a)+
  labs(x="variable",title = "Variables' importance score (METTL16, cut-off>0.6)")+
  theme_grey(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(size=18))
dev.off()



pdf("AUC_scat_Gen_06_16.pdf",width=7,height = 5) 
a<-as.integer(as.logical(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06)))
a[which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))]<-a[which(AUC_Gen_RF_16v_06==max(AUC_Gen_RF_16v_06))]-1
a<-a-7
a=a+8
ggplot(data=AUC_Gen_RF_16_df_06, aes(x=number, y=AUC))+
  
  geom_point(size=2,color=a)+
  labs(x="log(number of variables)",title = "AUC score over number of variables (METTL16, cut-off>0.6)")+
  theme(plot.title = element_text(size=15))
dev.off()
