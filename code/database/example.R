################# package #################


library(dplyr)
library(m6ALogisticModel)
library(ggplot2)
library(ROCR)
library(pROC)
library(e1071)

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

#To download package m6ALogisticModel, You need command: 
devtools::install_github("ZhenWei10/m6ALogisticModel")
#some package from bioconduct,you need download them by commmand like: 
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome")


################## Using data in SummarizedExperiment ##############

# you need your personal path 
SE <- readRDS("/media/NAS/kunqidir/DEseq/DEseq_06.rds")
 
# get methylation sets
# each row means a methylation set

rowRanges(SE)



#the information of relabtive enzyme (writer, earser)
# for example, to get the information about the "METTL3-consistent"


assay(SE[,grep("METTL3-consistent",colData(SE)$ID)])

# you can you use colData(SE)$ID to find which information you may need


####################### biological features generation ##################

library(m6ALogisticModel)

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

matureFE <- predictors_annot(se = SE,
                                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                   bsgnm = Hsapiens,
                                   fc = fitCons.UCSC.hg19,
                                   pc = phastCons100way.UCSC.hg19,
                                   struct_hybridize = Struc_hg19,
                                   feature_lst = Additional_features_hg19,
                                   hk_genes_list = HK_hg19_eids,
                                   genes_ambiguity_method = "average")

mcols(matureFE) # the each columns is a kind of bilogical features, and each row is a methylation site  

############################## sequence derived features ###############

# the sequence features from MethyRNA 

## get length of the sequences is 41 bp with the m6A motif in the center

seq_test <- as.character(DNAStringSet(Views(Hsapiens,rowRanges(SE)+20)))

## encoding 
source("Binary.R")

alter_chemicalNF(seq_test)


