# count_data shall be a matrix where rownames=proteins and colnames=runs
# annot_data shall be a data.frame with 3 columns: Run, Subject, Group

convert_to_MSstats <- function(count_data, annot_data) {
  count_data<-as.matrix(count_data)
  rownames(count_data)<-count_data[,1] # first column is protein name
  count_data<-count_data[,-1]
  rownames<-rownames(count_data)
  annot_data<-annot_data[,c("Run", "Subject", "Group")]
  # count_data <- prot_abundance
  # annot_data <- sample_annotation
  # 
  #Colnames of annotation file should be "Run", "Subject" and "Group"
  colnames(annot_data)[colnames(annot_data)=="Run"] <- "originalRUN"
  colnames(annot_data)[colnames(annot_data)=="Subject"] <- "SUBJECT_ORIGINAL"
  colnames(annot_data)[colnames(annot_data)=="Group"] <- "GROUP_ORIGINAL"
  
  # Melt spectral_count
  merged<-melt(count_data, value.name = "LogIntensities", varnames = c('Protein', 'originalRUN'))
  merged$LogIntensities <- suppressWarnings(as.numeric(paste(merged$LogIntensities)))
  
  # Merge with sample_annotation
  merged<-merge(merged, annot_data, by="originalRUN")
  
  # Set SUBJECT = SUBJECT_ORIGINAL
  merged[,"SUBJECT"] <- as.numeric(merged[, "SUBJECT_ORIGINAL"])
  
  # Set RUN = originalRUN
  merged[,"RUN"] <- as.numeric(merged[, "originalRUN"])
  
  # Set GROUP = GROUP_ORIGINAl
  merged[,"GROUP"] <- as.numeric(merged[, "GROUP_ORIGINAL"])
  
  # Create SUBJECT_NESTED
  merged[,"SUBJECT_NESTED"]<-paste(merged[,"SUBJECT"], merged[, "GROUP"], sep = ".")
  
  # Reorder columns to match example
  RunlevelData<-merged[,c("RUN", "Protein", "LogIntensities", "originalRUN", "GROUP", "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT")]
  
  # Remove samples with no condition specified
  RunlevelData<-na.omit(RunlevelData, cols="GROUP")
  ## make fake feature level data
  ProcessedData <- data.frame(PROTEIN = RunlevelData$Protein,
                              PEPTIDE = RunlevelData$Protein,
                              TRANSITION = RunlevelData$Protein,
                              FEATURE = RunlevelData$Protein,
                              LABEL = "H", 
                              GROUP_ORIGINAL = RunlevelData$GROUP_ORIGINAL,
                              SUBJECT_ORIGINAL = RunlevelData$SUBJECT_ORIGINAL,
                              RUN = RunlevelData$RUN,
                              GROUP = RunlevelData$GROUP,
                              SUBJECT = RunlevelData$SUBJECT,
                              INTENSITY = 2^RunlevelData$LogIntensities,
                              SUBJECT_NESTED = RunlevelData$SUBJECT_NESTED,
                              ABUNDANCE = RunlevelData$LogIntensities,
                              FRACTION = 1,
                              originalRUN = RunlevelData$originalRUN,
                              censored = FALSE,
                              SuggestToFilter = 0)
  
  QuantData <- list(ProcessedData = ProcessedData, 
                    RunlevelData = RunlevelData,
                    SummaryMethod = "TMP",
                    ModelQC = NULL,
                    PredictBySurvival = NULL)
  return(QuantData)
}
