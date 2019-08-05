#' Estimate the optimal size of training data for classification problem
#'
#' @param data Protein abundance data matrix. Rows are proteins and columns are samples.
#' @param group Group information for samples in data. The length of group should be equal to the number of columns in data.
#' @param n_sample number of different sample size to simulate. Default is 5 
#' @param sample_incr number of samples per condition to increase at each step. Default is 20 
#' @param protein_desc the fraction of proteins to reduce at each step. Proteins are ranked based on their mean abundance across all the samples. Default is 0.2. If protein_desc = 0.0, protein number will not be changed. 
#' @param iter number of times to repeat simulation experiments. Default is 10
#' @description For classification problem (such as disgnosys of disease), calculate the mean predictive accuray under different size of training data for future experiments of a Selected Reaction Monitoring (SRM), Data-Dependent Acquisition (DDA or shotgun), and Data-Independent Acquisition (DIA or SWATH-MS) experiment based on simulation.
#' @details The function fits intensity-based linear model on the input prelimiary data \emph{data} and uses variance components and mean abundance to simulate new training data with different sample size and protein number. Random forest model is fitted on simulated train data and used to predict the input preliminary data \emph{data}. The above procedure is repeated \emph{iter} times. Mean predictive accuracy and variance under different size of training data are reported. 
#' @return \emph{meanPA} is the mean predictive accuracy matrix under different size of training data.
#' @return \emph{varPA} is variance of predictive accuracy under different size of training data.
#' @author Ting Huang, Meena Choi, Olga Vitek.
#' @importFrom randomForest randomForest combine
#' 
#' @export
designSampleSizeClassificationCPTAC <- function(data, 
                                           group,
                                           n_sample = 5, 
                                           sample_incr = 20,
                                           protein_desc = 0.2,
                                           iter = 10,
                                           classifier = "randomForest") {
  
  message(" Preparing simulation...")

  ## Estimate the mean abundance and variance of each protein in each phenotype group
  parameters <- .eatimateBioVar(data,group)
 
  ## Prepare the parameters for simulation experiment
  mu <- parameters$mu
  sigma <- parameters$sigma
  promean <- parameters$promean
  valid_x <- as.data.frame(apply(data,1, function(x) .random.imp(x))) # impute missing values
  valid_y <- as.factor(group)
  
  ## Generate the vector of training sample size to simulate
  ngroup <- length(unique(group)) # Number of phenotype groups
  train_size <- seq.int(from = sample_incr, to = sample_incr * n_sample, length.out = n_sample)
  train_size <- train_size*ngroup
  message(" Size of training data to simulate: ", paste(train_size, collapse = ', '))
  
  ## Generate the vector of protein number to simulate
  nproteins <- nrow(mu)
  if(protein_desc == 0.0){ # no change of protein number
    protein_num <- nproteins
  } else{ # decrease the protein number by nproteins*protein_desc each time
    m_prot <- round(nproteins*protein_desc)
    protein_num <- seq.int(from = m_prot, to = nproteins, by = m_prot)
  }
  message(" Number of proteins to simulate: ", paste(protein_num, collapse = ', '))
  
  message(" Start to run the simulation...")
  
  PA <- list()
  for (i in 1:iter) { ## Number of iterations
    message("  Iteration: ", i)
    accur <- matrix(rep(0, times=length(train_size) * length(protein_num)), nrow=length(protein_num))
    
    ## simulate train data with different size
    for (n in seq_along(protein_num)) { 
      message("    Protein Number: ", protein_num[n])
      
      ## select proteins based on their mean abundance
      selectedPros<-order(promean,decreasing = TRUE)[1:protein_num[n]]
      mu_2 <- mu[selectedPros,]
      sigma_2 <-sigma[selectedPros,]
      
      for (m in seq_along(train_size)) { ## simulate samples in the training data
        ##Simulate training data 
        train <- .sampleSimulation(train_size[m], mu_2, sigma_2) 
        train_last <<- train
        x <- as.data.frame(train$X)
        colnames(x) <- rownames(mu_2)
        y <- as.factor(train$Y)
        ## Train random forest on training data
        rf <- randomForest::randomForest(x=x, y=y, mtry = )
        # fitControl <- trainControl(## 10-fold CV
        #   method = "repeatedcv",
        #   number = 10,
        #   ## repeated ten times
        #   repeats = 5)
        # rf <- train(x=x, y=y,
        #             method = "rf"
        #                  )
        ## Calculate predictive accuracy on validation data
        rf.pred <- predict(rf, valid_x) #Predict validation data
        accuracy <- sum(diag(table(rf.pred,valid_y))) / length(rf.pred)
        accur[n,m] <- accuracy
        
        # switch(classifier, 
        #        "randomForest"= {
        #           #Train random forest on training data
        #           rf <- randomForest::randomForest(x=x, y=y, mtry = )
        #           
        #           # rf<-train(y ~ ., 
        #           #           data=x,
        #           #           method="rf")
        #           ## Retired code: parallel computing for random forest
        #           #registerDoSNOW()
        #           #mcoptions <- list(set.seed=FALSE)
        #           #rf <- foreach(ntree=rep(10, 10), .combine=randomForest::combine, .multicombine=TRUE,
        #           #.packages='randomForest', .options.multicore=mcoptions) %dopar% {
        #           # randomForest(x=x,
        #           #y=y,
        #           #ntree=ntree)}
        #           
        #           ## Calculate predictive accuracy on validation data
        #           rf.pred <- predict(rf, valid_x)
        #           #Predict validation data
        #           accuracy <- sum(diag(table(rf.pred,valid_y))) / length(rf.pred)
        #           accur[n,m] <- accuracy
        #         },
        #        
        #         "xgBoost"= {
        #           message("Classifier not implemented")
        #           return()
        #         },
        #         {
        #           message("No classifier selected")
        #           return()
        #         }
        #          )
        
      }
    }
    PA[[i]] <- accur
  }
  
  ## Calculate the mean accuracy and variance
  meanPA <- matrix(rep(0, times=length(train_size) * length(protein_num)), nrow=length(protein_num))
  varPA <- matrix(rep(0, times=length(train_size) * length(protein_num)), nrow=length(protein_num))
  
  for (n in seq_along(protein_num)) {
    for (m in seq_along(train_size)) {
      temp <- NULL
      for (i in 1:iter) {
        temp <- c(temp, PA[[i]][n, m])
      }
      meanPA[n, m] <- mean(temp)
      varPA[n, m] <- var(temp)
    }
  }
  rownames(meanPA) <- paste0("prot", protein_num)
  colnames(meanPA) <- paste0("tra", train_size)
  rownames(varPA) <- paste0("prot", protein_num)
  colnames(varPA) <- paste0("tra", train_size)
  message(" Simulation completed.")
  
  return(list(meanPA = meanPA, 
              varPA = varPA))
}


#' Simulate extended datasets for sample size estimation
#'
#' @param m number of samples to simulate
#' @param mu a matrix of mean abundance in each phenotype group and each protein
#' @param sigma a matrix of variance in each phenotype group and each protein
#' @return \emph{X} A simulated matrix with required sample size
#' @return \emph{Y} Group information corresponding with \emph{X}
#' @keywords internal
.sampleSimulation <- function(m, mu, sigma) {
  
  nproteins <- nrow(mu)
  ngroup <- ncol(mu)
  ## Determine the size of each phenotype group
  samplesize <- .determineSampleSizeinEachGroup(m, ngroup)
  
  ## Simulate the data matrix
  sim_matrix <- matrix(rep(0, nproteins * m), ncol=m)
  for (i in 1:nproteins) {
    abun <- NULL
    for (j in 1:ngroup) {
      abun <- c(abun, rnorm(samplesize[j], mu[i, j], sigma[i, j]))
    }
    sim_matrix[i, ] <- abun
  }
  sim_matrix <- t(sim_matrix)
  colnames(sim_matrix) <- rownames(mu)
  #Simulate the phenotype information
  group <- rep(c(1:ngroup), times=samplesize)
  
  index <- sample(length(group), length(group))
  sim_matrix <- sim_matrix[index, ]
  group <- group[index]
  
  return(list(X=sim_matrix,
              Y=as.factor(group)))
}

#' Determine the size of each phenotype group
#'
#' @param m sample size
#' @param ngroup number of phenotype groups
#' @return vector of sample size in each group
#' @keywords internal
.determineSampleSizeinEachGroup <- function(m, ngroup) {
  samplesize <- vector(mode="numeric", length=ngroup)
  counter <- 0
  
  while (counter < m) {
    for (i in 1:ngroup) {
      if (counter < m) {
        counter <- counter + 1
        samplesize[i] <- samplesize[i] + 1
      }
    }
  }
  return(samplesize)
}

#' Estimate the mean abundance and variance of each protein in each phenotype group.
#'
#' @param data Protein abundance data matrix.
#' @return \emph{mu} is the mean abundance matrix of each protein in each phenotype group; 
#' @return \emph{sigma} is the sd matrix of each protein in each phenotype group;#' @keywords internal
.eatimateBioVar <-  function(data,group) {
  biovar<-NULL
  meanabun<-NULL
  ngroup <- length(unique(group))
  
  GroupVar <- matrix(rep(NA, nrow(data) * ngroup), ncol = ngroup)
  GroupMean <- matrix(rep(NA, nrow(data) * ngroup), ncol = ngroup)
  SampleMean <- NULL # mean across all the samples
  for (i in 1:nrow(data)) {
    sub<- NULL
    sub$ABUNDANCE <- data[i,]
    sub$GROUP <- factor(group)
    sub<-do.call(cbind.data.frame, sub)
    sub <- sub[!is.na(sub$ABUNDANCE),]
    
    df.full <- lm(ABUNDANCE ~ GROUP , data = sub)
    var<-anova(df.full)[[3]][2]
    
    abun<-summary(df.full)$coefficients[,1]
    abun[-1]<-abun[1]+abun[-1]
    GroupVar[i,]<-rep(sqrt(var),times=length(abun))
    GroupMean[i,] <-abun
    SampleMean <- c(SampleMean, mean(sub$ABUNDANCE, na.rm = T))
  }
  
  rownames(GroupVar) <- rownames(data)
  rownames(GroupMean) <- rownames(data)
  names(SampleMean) <- rownames(data)
  return(list(mu=GroupMean,sigma=GroupVar, promean = SampleMean))
}

#' For each protein, impute the missing values based on the observed values
#'
#' @param data protein abundance data for one protein.
#' @return Imputed protein abundance data
#' @keywords internal
.random.imp <- function (data){
  missing <- is.na(data) # count missing values
  n.missing <- sum(missing)
  data.obs <- data[!missing] # keep the observed values
  imputed <- data
  # impute missing values by randomly selecting observed values
  imputed[missing] <- sample (data.obs, n.missing, replace=TRUE)
  return (imputed)
} 

#############################################
## designSampleSizeClassificationPlots
#############################################
#' Visualization for sample size calculation in classification problem
#' @param data output from function \code{\link{designSampleSizeClassification}}
#' @description To illustrate the mean classification accuracy under different protein number and sample size. The input is the result from function \code{\link{designSampleSizeClassification}}.
#' @return Plot for sample size estimation. x-axis : sample size, y-axis: mean predictive accuracy. Color: different protein number.
#' @details Data in the example is based on the results of sample size calculation in classification problem from function \code{\link{designSampleSizeClassification}}
#' @author Ting Huang, Meena Choi, Olga Vitek. 
#' @importFrom reshape2 melt
#' @export
designSampleSizeClassificationPlots <- function(data) {
  
  ## ggplot needs a long format dataframe
  ## get the mean accuracy
  meandata <- as.data.frame(data$meanPA)
  meandata$Protein_number <- rownames(meandata)
  meandata <- reshape2::melt(meandata, id.vars = "Protein_number", variable.name = "Train_size", value.name = "mean")
  
  ## get the variance
  vardata <- as.data.frame(data$varPA)
  vardata$Protein_number <- rownames(vardata)
  vardata <- reshape2::melt(vardata, id.vars = "Protein_number", variable.name = "Train_size", value.name = "var")
  
  ## perform the join
  plotdata <- merge(meandata, vardata, all=TRUE)
  # get standard deviation column
  plotdata$sd <- sqrt(plotdata$var)
  # make sure train size is numeric
  plotdata$Train_size <- gsub("tra", "", plotdata$Train_size)
  plotdata$Train_size <- as.numeric(as.character(plotdata$Train_size))  
  # make sure Protein_number is ordered factor
  plotdata$Protein_number <- gsub("prot", "", plotdata$Protein_number)
  plotdata$Protein_number <- factor(plotdata$Protein_number, levels = sort(as.numeric(unique(plotdata$Protein_number))))
  
  ## make the plot
  p <- ggplot(data = plotdata, aes(x= Train_size, y= mean, group = Protein_number, colour = Protein_number)) +
    geom_point() +
    geom_line() +
    labs(title="Sample size estimation", x="Sample size", y='Mean accuracy') +
    guides(color=guide_legend(title="Protein Number"))
  
  print(p)
}