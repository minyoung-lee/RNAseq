rm(list=ls())
setwd("E:/miRNA_NCI60/work_matlab")
library("limma")

filenameX='LIMMA_inputX.txt'
dataX <- as.matrix(read.table(filenameX, header=F, row.names=NULL, sep="\t", check.names=F))
filenameY='LIMMA_inputY.txt'
dataY=as.matrix(read.table(filenameY, header=F, row.names=NULL, sep="\t", check.names=F))


# assume dataX has columnes for samples and features as row
#save_name="EM"
p_cutoff=0.05
FC_cutoff=1.5

  # build the design matrix
  control <- rep(0, length(dataY))
  disease <- rep(0, length(dataY))
  control[dataY == "control"] <- 1
  disease[dataY == "disease"] <- 1
  design <- cbind(control, disease)
  colnames(design) <- c("control", "disease")
  
  # fit the model
  fit <- lmFit(dataX, design)
  cont.matrix <- makeContrasts(diseasevscontrol="disease-control", levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  result <- topTable(fit2, number="all", adjust.method = "BH", sort.by="none")
  #write.table(result, paste0("limma_results_",save_name,".txt"), sep="\t", col.names=NA)
write.table(result, "limma_results.txt", sep="\t", col.names=NA)



