simulatedClassificationCrossClassifiers <- function(data,
                                                    parameters, 
                                                    n_sample = 5, 
                                                    n_protein = 20, 
                                                    iter = 10,
                                                    classifiers = c("rf", "nnet")
) {
  
  ## Save a copy of the parameters used
  paramlog <- c("n_sample"=n_sample, "sample_incr" = sample_incr, "protein_desc" = protein_desc, "iter" = iter, "classifier" = classifier)
  
  ## Prepare the parameters for simulation experiment
  mu <- parameters$mu
  sigma <- parameters$sigma
  promean <- parameters$promean
  valid_x <- as.data.frame(parameters$X)
  debug_valid_x<<-valid_x
  valid_y <- as.factor(parameters$Y)
  debug_valid_y<<-valid_y
  
  message(" Datatset parameters: ", paste(n_sample, n_protein, collapse = ', '))
  
  message(" Classifiers: ", classifiers)
  
  message(" Start to run the simulation...")
  
  PA <- list()
  models <- list()
  for (i in 1:iter) { ## Number of iterations
    message("  Iteration: ", i)
    accur <- matrix(rep(0, times=length(classifier)))
    models[[i]]<-list()
    ## simulate train data with different size
    for (classifier in classifiers) { ## simulate samples in the training data
      ##Simulate training data 
      selectedPros<-order(promean,decreasing = TRUE)[1:n_protein]
      mu_2 <- mu[selectedPros,]
      sigma_2 <-sigma[selectedPros,]
      
      train <- .sampleSimulation(n_sample, mu_2, sigma_2) 
      x <- as.data.frame(train$X)
      colnames(x) <- rownames(mu_2)
      y <- as.factor(train$Y)
      
      start_time <- Sys.time()
      
      #Train random forest on training data
      
      cl <- makePSOCKcluster(5)
      registerDoParallel(cl)
      switch(classifier,
             "rf"= {tunegrid=data.frame(mtry=2)},
             "nnet"={tunegrid=data.frame(size=5, decay=0.1)},
             "svmLinear"={tunegrid=data.frame(C=1)},
             "pls"={tunegrid=data.frame(ncomp=2)},
             "naive_bayes"={tunegrid=data.frame(laplace=0, usekernel=FALSE, adjust=1)}
      )
      model <- .train_model(x, y, classifier, tunegrid)
      # if (use_caret==TRUE) {stopCluster(cl)}
      stopCluster(cl)
      
      end_time <- Sys.time()
      print(end_time - start_time)
      
      ## Calculate predictive accuracy on validation data
      model.pred <- predict(model, valid_x) #Predict validation data
      accuracy <- sum(diag(table(model.pred,valid_y))) / length(model.pred)
      accur[classifier] <- accuracy
      models[[i]][[classifier]]<-model
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
  
  ##
  processout <- rbind(processout, 
                      "The number of sample and proteins is estimated - okay")
  write.table(processout, file=finalfile, row.names=FALSE)
  
  return(list(meanPA = meanPA, 
              varPA = varPA,
              models = models,
              paramlog = paramlog))
}