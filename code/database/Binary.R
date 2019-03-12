
require(dplyr)
motifSplit <- function(data){
  motifSplit <- list(1)
  for (i in 1:length(data)){
    motifstr <- as.vector(strsplit(data[i],"",fixed = TRUE))
    motifSplit <- c(motifSplit,motifstr)
  }
  motifSplit[[1]] <- NULL
  return(motifSplit)
} 

#motiftest <- motifSplit(motif)

binary_encoding <- function(Data){
  A <- c(1,0,0,0)%>%matrix(nrow=1)
  C <- c(0,1,0,0)%>%matrix(nrow=1)
  G <- c(0,0,1,0)%>%matrix(nrow=1)
  T <- c(0,0,0,1)%>%matrix(nrow=1)
  binary_matrix <- matrix(NA,nrow = 1, ncol = 1)
  for(i in 1:length(Data)){
    if (Data[i]=="A"){
      binary_matrix <- cbind(binary_matrix,A)
    }else if (Data[i]=="C"){
      binary_matrix <- cbind(binary_matrix,C)
    }else if (Data[i]=="G"){
      binary_matrix <- cbind(binary_matrix,G)
    }else if (Data[i]=="T"){
      binary_matrix <- cbind(binary_matrix,T)
    }
  }
  binary_matrix <- binary_matrix[,-1]
  return(binary_matrix)
}

motifBinary <- function(data){
  motiftest <- motifSplit(data)
  motifBinary <- matrix(NA,nrow = 1, ncol = 4*length(motiftest[[1]]))
  for (i in 1:length(motiftest)){
    motifstr <- binary_encoding(unlist(motiftest[[i]]))%>%matrix(nrow = 1)
    motifBinary <- rbind(motifBinary,motifstr)
  }
  motifBinary <- motifBinary[-1,]
  return(motifBinary)
}

alter_chemical_encoding <- function(Data){
  A <- c(1,1,1)%>%matrix(nrow=1)
  C <- c(0,1,0)%>%matrix(nrow=1)
  G <- c(1,0,0)%>%matrix(nrow=1)
  T <- c(0,0,1)%>%matrix(nrow=1)
  binary_matrix <- matrix(NA,nrow = 1, ncol = 1)
  for(i in 1:length(Data)){
    if (Data[i]=="A"){
      binary_matrix <- cbind(binary_matrix,A)
    }else if (Data[i]=="C"){
      binary_matrix <- cbind(binary_matrix,C)
    }else if (Data[i]=="G"){
      binary_matrix <- cbind(binary_matrix,G)
    }else if (Data[i]=="T"){
      binary_matrix <- cbind(binary_matrix,T)
    }
  }
  binary_matrix <- binary_matrix[,-1]
  return(binary_matrix)
}

alter_chemical <- function(data){
  motiftest <- motifSplit(data)
  motifBinary <- matrix(NA,nrow = 1, ncol = 3*length(motiftest[[1]]))
  for (i in 1:length(motiftest)){
    motifstr <- alter_chemical_encoding(unlist(motiftest[[i]]))%>%matrix(nrow = 1)
    motifBinary <- rbind(motifBinary,motifstr)
  }
  motifBinary <- motifBinary[-1,]
  return(motifBinary)
}

alter_chemical_NF <- function(Data){
  A <- c(1,1,1)%>%matrix(nrow=1)
  C <- c(0,1,0)%>%matrix(nrow=1)
  G <- c(1,0,0)%>%matrix(nrow=1)
  T <- c(0,0,1)%>%matrix(nrow=1)
  N <- c(0,0,0,0)%>%matrix(nrow=1)
  binary_matrix <- matrix(NA,nrow = 1, ncol = 1)
  for(i in 1:length(Data)){
    if (Data[i]=="A"){
      binary_matrix <- cbind(binary_matrix,A)
      TTAppear <- length(which(Data[1:i]=="A"))
      Nucleotide_frequency <- (TTAppear/i) 
      binary_matrix <- cbind(binary_matrix,Nucleotide_frequency)
    }else if (Data[i]=="C"){
      binary_matrix <- cbind(binary_matrix,C)
      TTAppear <- length(which(Data[1:i]=="C"))
      Nucleotide_frequency <- (TTAppear/i) 
      binary_matrix <- cbind(binary_matrix,Nucleotide_frequency)
    }else if (Data[i]=="G"){
      binary_matrix <- cbind(binary_matrix,G)
      TTAppear <- length(which(Data[1:i]=="G"))
      Nucleotide_frequency <- (TTAppear/i) 
      binary_matrix <- cbind(binary_matrix,Nucleotide_frequency)
    }else if (Data[i]=="T"){
      binary_matrix <- cbind(binary_matrix,T)
      TTAppear <- length(which(Data[1:i]=="T"))
      Nucleotide_frequency <- (TTAppear/i) 
      binary_matrix <- cbind(binary_matrix,Nucleotide_frequency)
    }else{
      binary_matrix <- cbind(binary_matrix,N)
    }
  }
  binary_matrix <- binary_matrix[,-1]
  return(binary_matrix)
}

alter_chemicalNF <- function(data){
  motiftest <- motifSplit(data)
  motifBinary <- matrix(NA,nrow = 1, ncol = 4*length(motiftest[[1]]))
  for (i in 1:length(motiftest)){
    motifstr <- alter_chemical_NF(unlist(motiftest[[i]]))%>%matrix(nrow = 1)
    motifBinary <- rbind(motifBinary,motifstr)
  }
  motifBinary <- motifBinary[-1,]
  return(motifBinary)
}

Accuracy <- function(data){
  Accury <- (data[1,1]+data[2,2])/(data[1,2]+data[2,2]+data[1,1]+data[2,1])
  return(Accury)
} 

Specificity <- function(data){
  Specificity <- (data[2,2])/(data[2,1]+data[2,2])
  return(Specificity)
}

Sensitiveity <- function(data){
  Sensitiveity <- (data[1,1])/(data[1,2]+data[1,1])
  return(Sensitiveity)
}

MCC <- function(data){
  MCC <- (data[1,1]*data[2,2]-data[1,2]*data[2,1])/
    (sqrt(as.numeric(data[1,1]+data[1,2])*(data[1,1]+data[2,1])*(data[2,2]+data[2,1])*(data[2,2]+data[1,2])))
  return(MCC)
} 

F1 <-  function(data){
  F1 <- (2*data[1,1])/(data[1,2]+2*data[1,1]+data[2,1])
  return(F1)
} 

ML_ALL <- function(data){
  ML_ALL <- matrix(c(Accuracy(data), Sensitiveity(data),Specificity(data), MCC(data), F1(data),
                     "Accuracy","Sensitiveity","Specificity","MCC","F1"),ncol = 2)
  return(ML_ALL)
}
