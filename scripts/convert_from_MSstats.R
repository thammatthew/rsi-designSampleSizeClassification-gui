library(MSstats)
library(MSstatsBioData)
library(reshape2)

data_processed <- dataProcess(SRM_crc_training)

prot_abundance <- as.matrix(dcast(data_processed[["RunlevelData"]], Protein~SUBJECT_ORIGINAL, value.var = "LogIntensities"))
rownames(prot_abundance)<-prot_abundance[,1] # first column is protein name
prot_abundance<-prot_abundance[,-1]
sample_annotation <- data_processed[["RunlevelData"]][,c("SUBJECT_ORIGINAL", "SUBJECT_ORIGINAL","GROUP_ORIGINAL")][!duplicated(data_processed[["RunlevelData"]][,c("SUBJECT_ORIGINAL", "SUBJECT_ORIGINAL","GROUP_ORIGINAL")]),]
colnames(sample_annotation) <- c("Run", "Subject", "Group")

