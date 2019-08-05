function(session, input, output) {
  # Set maximum size of uploaded files to 300mb
  options(shiny.maxRequestSize=300*1024^2)
  # Create reactiveValue object
  values <- reactiveValues()
  # Reactive expression to check status of designSampleSizeClassification
  output$sim <- reactive({!is.null(values$result)})
  outputOptions(output, "sim", suspendWhenHidden = FALSE)
  
  # Debug functions
  # output$res<-renderText({
  #   paste("You've selected:", input$tabs)
  # })
  observeEvent(input$debug_save, {
    try ({
      message("Starting to output environment")
      
      # T1+2 Data
      debug_quant_data <<- values$quant_data
      debug_prot_abundance <<- values$prot_abundance
      debug_sample_annotation <<- values$sample_annotation
      debug_pca <<- values$pca
      debug_parameters <<- values$parameters
      debug_s_summary <<- values$s_summary
      
      # T3+4 Data
      debug_simulated_raw <<- values$simulated_raw
      debug_valid_y<<-values$valid_y
      debug_valid_x<<-values$valid_x
      debug_chosen_dataset<<-values$chosen_dataset
      debug_result <<- values$result
      debug_meandata<<- meandata
      debug_chosen_model <<- values$chosen_model
      debug_chosen_model.pred <<- values$chosen_model.pred
      debug_confusionMatrix <<- values$confusionMatrix
      debug_chosen_model.metrics <<- values$chosen_model.metrics
      debug_gg_meanPA <<- values$gg_meanPA
      message("Environment output complete")
      showNotification("Environment Saved", duration=10, closeButton = TRUE, type="message")
      
    })
  })
  
  # T1 Sidebar Data Import
  ## L1 Options Menu
  output$options<-renderUI({
    if (is.null(input$data_format)) {
      tags$b("Please select a data format")
    } else {
      switch(input$data_format,
             "standard"= tagList(
               selectInput("separator", "Separator", choices = list("Commas"=",", "Tabs"="\t"))
             ),
             "MSstats"=tagList(
               tags$b("placeholder")
             ),
             "Examples"=tagList(
               tags$b("placeholder")
             ),
             tags$b("No available options for selected data format.")
            )
    }
  })  
  outputOptions(output, "options", suspendWhenHidden = FALSE)

  ## L1 Advanced Options Menu 
  output$advanced_options<-renderUI({
    if (is.null(input$data_format)) {
      tags$b("Please select a data format")
    } else {
      switch(input$data_format,
             "standard"= tagList(
               tags$b("Preprocessing"),
               checkboxInput("cptac_log2", label="Log2 Transformation", value=FALSE),
               checkboxInput("standard_quantnorm", label="Quantile Normalisation", value=FALSE),
               sliderInput("missingthreshold", label = tags$b("Threshold for Missing Values"), min = 0, max = 100, value = 10)
             ),
             "MSstats"=tagList(
               tags$b("placeholder")
             ),
             "Examples"=tagList(
               tags$b("placeholder")
             ),
             tags$b("No available options for selected data format.")
      )
    }
  })
  outputOptions(output, "advanced_options", suspendWhenHidden = FALSE)
  
  ## L1 Dataset Upload and Import Menu | Contingent on Options & Advanced Options loading
  output$select_files<-renderUI({
    if (is.null(input$data_format)) {
      tags$b("Please select a data format")
    } else {
      switch(input$data_format,
             "standard" = tagList(
                 fileInput("standard_count", "Select Protein Abundance File",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv", "text/tab-separated-values", ".tsv")),
                 fileInput("standard_annot", "Select Sample Annotation File",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv", "text/tab-separated-values", ".tsv")),
                 actionButton("import_data", "Import and Preprocess Data", icon=icon("file-import"))
               ),
             "MSstats" = tagList(
               fileInput("msstats_object", "Select MSstats-processed .rda Object",
                         multiple = FALSE,
                         accept = c(".rda", ".rdata")),
               actionButton("import_data", "Import and Preprocess Data", icon=icon("file-import"))
             ),
             "Examples" = tagList(
               selectInput("example_dataset", "Select an Example Dataset", choices = data(package = "MSstatsBioData")$result[, "Item"]),
               actionButton("import_data", "Import and Preprocess Data", icon=icon("file-import"))
             ),
             "Debug" = tagList(
               fileInput("environment_image", "Select Environment Image",
                         multiple = FALSE),
               actionButton("import_data", "Import and Preprocess Data", icon=icon("file-import"))
             ),
             tags$b("No import method for selected data format.")
      )
    }
  })
  outputOptions(output, "select_files", suspendWhenHidden = FALSE)
  
  # Actual Data import function
  observeEvent(input$import_data, {
    if (is.null(input$data_format)) {
      showNotification("No imported dataset found", duration=10, closeButton = TRUE, type="message")
    } else if (input$data_format=="Debug") {
      load(input$environment_image$datapath)
      
      ### Legacy code from previous variant
      # values$data<-debug_data
      # values$quant_data<-debug_data_processed
      # values$prot_abundance<-debug_prot_abundance
      # values$sample_annotation<-debug_sample_annotation
      # 
      # tprot <- debug_tprot2
      # tprot_pca1 <- prcomp(tprot, center = TRUE,scale. = TRUE)
      # values$pca<-tprot_pca1
      # showNotification("PCA complete")
      # 
      # values$result<-debug_result
      # values$meandata<-debug_meandata
      # values$gg_meanPA<-melt(debug_result$meanPA)
      # # values$gg_meanPA$Var1<-as.character(values$gg_meanPA$Var1)
      # # values$gg_meanPA$Var2<-as.character(values$gg_meanPA$Var2)
      values$dataset_name<-input$environment_image$name
      values$is_imported<-TRUE
    } else {
      switch(input$data_format,
             "standard" = {
               tryCatch(
                 {
                   # Load data from user-input filepath
                   values$count_data <- read.csv(input$standard_count$datapath)
                   values$annot_data <- read.csv(input$standard_annot$datapath)
                   values$quant_data <- convert_to_MSstats(values$count_data, values$annot_data)
                   values$dataset_name <- input$standard_count$name
                   values$is_imported <- TRUE
                   showNotification("Data import complete", duration=10, closeButton = TRUE, type="message")
                 },
                 error = function(e) {
                   showNotification("Error", duration=10, closeButton = TRUE, type="message")
                 }
               )
             },
             
             "MSstats" = {
               tryCatch(
                 {
                   values$quant_data<-load(input$msstats_object$datapath)
                   values$dataset_name<-input$msstats_object$name
                   values$is_imported<-TRUE
                   showNotification("Data import complete", duration=10, closeButton = TRUE, type="message")
                 },
                 error = function(e) {
                   showNotification("Error", duration=10, closeButton = TRUE, type="message")
                 }
               )
             },
             "Examples" = {
               tryCatch(
                 {
                   values$data <- get(input$example_dataset)
                   values$quant_data <- dataProcess(values$data)
                   values$dataset_name<-as.character(input$example_dataset)
                   values$is_imported<-TRUE
                   showNotification("Data import complete", duration=10, closeButton = TRUE, type="message")
                 },
                 error = function(e) {
                   showNotification("Error", duration=10, closeButton = TRUE, type="message")
                 }
               )
             }
      )
      prot_abundance <- as.matrix(dcast(values$quant_data[["RunlevelData"]], Protein~SUBJECT_ORIGINAL, value.var = "LogIntensities"))
      rownames(prot_abundance)<-prot_abundance[,1] # first column is protein name
      values$prot_abundance<-prot_abundance[,-1]
      values$sample_annotation <- values$quant_data[["RunlevelData"]][,c("SUBJECT_ORIGINAL","GROUP_ORIGINAL")]
      summary.s <- matrix(NA,ncol=nlevels(values$quant_data[["RunlevelData"]]$GROUP_ORIGINAL), nrow=3)
      ## # of MS runs
      temp <- unique(values$quant_data[["RunlevelData"]][, c("GROUP_ORIGINAL", "RUN")])
      temp1 <- xtabs(~GROUP_ORIGINAL, data=temp)
      summary.s[1,] <- temp1
      ## # of biological replicates
      temp <- unique(values$quant_data[["RunlevelData"]][, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL")])
      temp1 <- xtabs(~GROUP_ORIGINAL, data=temp)
      summary.s[2,] <- temp1
      ## # of technical replicates
      c.tech <- round(summary.s[1,] / (summary.s[2,] * length(unique(values$quant_data[["ProcessedData"]]$FRACTION))))
      ##summary.s[3,] <- ifelse(c.tech==1,0,c.tech)
      summary.s[3,] <- c.tech
      colnames(summary.s) <- unique(values$quant_data[["RunlevelData"]]$GROUP_ORIGINAL)
      rownames(summary.s) <- c("# of MS runs","# of Biological Replicates", "# of Technical Replicates")
      values$summary <- summary.s
      showNotification("Data processing complete", duration=10, closeButton = TRUE, type="message")
      
      tprot <- t(values$prot_abundance)
      class(tprot) <- "numeric"
      tprot[is.na(tprot)] <- 0

      remove_cols<-nearZeroVar(tprot, names=TRUE, freqCut=19, uniqueCut=10)
      keep_cols<-colnames(tprot)
      tprot<-tprot[,setdiff(keep_cols,remove_cols)]
      
      tprot_pca1 <- prcomp(tprot, center = TRUE,scale. = TRUE)
      values$pca<-tprot_pca1
      showNotification("PCA complete")
      
      values$parameters<-buildSimulatedDistribution(values$quant_data)
      showNotification("Analysis of variance complete")
    }
  })
  
  # T2 Explore Dataset Tab
  ## Main Content
  output$explore_data_content<-renderUI({
      if (is.null(values$is_imported)) {
        tagList(
          h1("Explore Dataset"),
          tags$b("Please use the 'Import Data' menu to import a proteome dataset.")
        )
      } else {
        n_prot <- nrow(unique(values$quant_data[["RunlevelData"]]["Protein"]))
        n_group <- nrow(unique(values$quant_data[["RunlevelData"]]["GROUP_ORIGINAL"]))
        n_biorep <- nrow(unique(values$quant_data[["RunlevelData"]]["SUBJECT"]))
        n_techrep <- floor(nrow(unique(values$quant_data[["RunlevelData"]]["originalRUN"]))/nrow(unique(values$quant_data[["RunlevelData"]]["SUBJECT"])))
        tagList(
          h1(paste("Explore",values$dataset_name)),
          fluidRow(
            valueBox("Proteins Quantified", value = n_prot, icon=icon("dna"), color="purple", width=6),
            valueBox("Groups", value = n_group, icon=icon("layer-group"), color="green", width=6)
          ),
          fluidRow(
            valueBox("Biological Replicate(s)", value = n_biorep, icon=icon("users"), color="aqua", width=6),
            valueBox("Technical Replicates(s)", value = n_techrep, icon=icon("copy"), color="light-blue", width=6)
            
          ),
          fluidRow(
            tabBox(id="tabset1", width=12,
                   tabPanel("Raw Data",
                            div(style = 'overflow-x:scroll;',
                                DT::dataTableOutput("sample_annotation_table"),
                                DT::dataTableOutput("prot_abundance_table")
                            )
                   ),
                   tabPanel("Summary",
                            div(style = 'overflow-x:scroll;',
                                DT::dataTableOutput("summary_table")
                            )
                   ),
                   tabPanel("PCA",
                            div(style = 'overflow-x:scroll;',
                                plotOutput("pcabiplot"),
                                plotOutput("pcascreeplot")
                            )
                   ),
                   tabPanel("QC Box Plots",
                            div(style = 'overflow-x:scroll;',
                                plotlyOutput("global_boxplot")
                            )
                   ),
                   tabPanel("Mean-variance Plots",
                            div(style = 'overflow-x:scroll;',
                                plotlyOutput("mean_variance_plot")))
            )
          )
        )
      }
  })
  
  #L1 Explore Dataset / Data Tables
  output$prot_abundance_table = DT::renderDataTable({
    values$prot_abundance
  })
  output$sample_annotation_table = DT::renderDataTable({
    values$sample_annotation
  })
  output$summary_table = DT::renderDataTable({
    values$summary
  })
  
  # L1 Explore Dataset / Plots
  ## PCA
  output$pcabiplot <- renderPlot({
    ggbiplot(values$pca, ellipse=TRUE, var.axes=FALSE, groups=values$sample_annotation[!duplicated(values$sample_annotation),][,2])
    # plot_ly(values$pca, x=Comp.1, y=Comp.2, text=rownames(values$pca), mode="markers", marker=list(size=11))
  })
  
  output$pcascreeplot <- renderPlot({
    ggscreeplot(values$pca)
  })
  
  ## Global QC Box Plot
  output$global_boxplot <- renderPlotly({
    data <- values$quant_data[["RunlevelData"]]
    ylim1 = boxplot.stats(data$y)$stats[c(1,5)]
    
    plot_ly(data, y=~LogIntensities, x=~originalRUN, color=~GROUP_ORIGINAL, type = "box") %>%
      layout(xaxis=list(title="Biological Replicate"), yaxis=list(range = ylim1, title="Log Intensities"))
    

  })
  
  ## Mean-variance Plot
  output$mean_variance_plot <- renderPlotly({
    data <- as.data.frame(cbind(mean = values$parameters$mu[,1], sd = values$parameters$sigma[,1]))
    plot.lowess <- lowess(data)
    plot.lowess <- data.frame(x = plot.lowess$x, y = plot.lowess$y)
    
    ggplot(data = data, aes(mean, sd)) +
      stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
      scale_fill_continuous(low = "white", high = "#0072B2")+
      geom_point(alpha = 0.02, shape = 20)+
      # geom_vline(xintercept=sim_quan[2:5], color="gray") +
      geom_line(data = plot.lowess, aes(x, y), color="orange", size = 1) +
      scale_y_continuous(expand = c(0,0), limits = c(0,3.01)) + 
      scale_x_continuous(expand = c(0,0), limits = c(0,30.5)) +
      theme_bw()+ 
      theme(legend.position="none", axis.text=element_text(size=15))+ 
      labs(x = "Mean log2 intensity")+
      labs(y = "Standard deviation")
    
  })
  
  # T3 Explore Simulated Data Tab
  output$explore_simulated_content<-renderUI({
    if (is.null(values$is_imported)) {
      tags$b("Please use the 'Import Data' menu to import a proteome dataset.")
    } else {
      tagList(
        fluidRow(
          box(title = paste("Selected Dataset:",values$selected_dataset_name), width=9,
              htmlOutput("inspect_simulated")
          ),
          box(width=3, title = "Simulated Datasets",
            plotlyOutput("simulated_grid"),
            br(),
            h3("Parameters", class="custom-box-title"),
            numericInput("n_sample", label="Number of Different Sample Sizes to Simulate", value=5),
            numericInput("sample_incr", label="Step Size Between Simulated Sample Sizes", value=20),
            numericInput("n_protein", label="Number of Different Protein Counts to Simulate", value=5)
          )
       )
      )
    }
  })
  outputOptions(output, "explore_simulated_content", suspendWhenHidden = FALSE)
  
  #L1 Explore Simulated Data / Process and generate tabset for specified simulated dataset
  output$inspect_simulated<-renderUI({
    values$selected_dataset <- event_data("plotly_click", source = "simulated_grid_plot")
    if (length(values$selected_dataset)) {
      protein_desc <- 1/input$n_protein
      
      mu <- values$parameters$mu
      sigma <- values$parameters$sigma
      promean <- values$parameters$promean

      # Use mean and variance components to generate dataset on the fly
      ## No. of proteins is vars[1], No. of samples is vars[2]
      vars <- c(values$selected_dataset[["x"]], values$selected_dataset[["y"]])
      values$selected_dataset_name<-paste(vars[1]," Proteins, ",vars[2]," BioReplicates",sep="")
      
      selectedProt<-order(promean, decreasing=TRUE)[1:vars[1]]
      mu_selected <- mu[selectedProt,]
      sigma_selected <- sigma[selectedProt,]
      simulated_data <- .sampleSimulation(as.integer(vars[2]), mu_selected, sigma_selected)

      values$s_prot_abundance <- t(simulated_data[["X"]])
      values$s_sample_annotation <- simulated_data[["Y"]]
      summary.s <- matrix(NA,ncol=nlevels(values$s_sample_annotation), nrow=1)
      
      merged<-melt(simulated_data[["X"]], value.name = "LogIntensities", varnames = c('originalRUN', 'Protein'))
      merged$LogIntensities <- suppressWarnings(as.numeric(paste(merged$LogIntensities)))
      values$s_groups<-data.frame("originalRUN"=1:length(simulated_data[["Y"]]), "Group"=simulated_data[["Y"]])
      merged<-merge(merged, values$s_groups, by="originalRUN")
      values$s_quant_data<-merged
      
      ## # of biological replicates
      temp <- unique(merged[, c("Group", "originalRUN")])
      temp1 <- xtabs(~Group, data=temp)
      summary.s[1,] <- temp1
      
      colnames(summary.s) <- unique(values$quant_data[["RunlevelData"]]$GROUP_ORIGINAL)
      rownames(summary.s) <- c("# of Biological Replicates")
      values$s_summary <- summary.s
      showNotification("Data processing complete", duration=10, closeButton = TRUE, type="message")
      
      tprot <- t(values$s_prot_abundance)
      class(tprot) <- "numeric"
      
      remove_cols<-nearZeroVar(tprot, names=TRUE, freqCut=19, uniqueCut=10)
      keep_cols<-colnames(tprot)
      tprot<-tprot[,setdiff(keep_cols,remove_cols)]
      
      tprot_pca1 <- prcomp(tprot, center = TRUE,scale. = TRUE)
      values$s_pca<-tprot_pca1
      showNotification("PCA complete")
      
      tabBox(id="simulated_data_wrapper", width="100%",
             tabPanel("Raw Data",
                      div(style = 'overflow-x:scroll;',
                          DT::dataTableOutput("s_sample_annotation_table"),
                          DT::dataTableOutput("s_prot_abundance_table")
                      )
             ),
             tabPanel("Summary",
                      div(style = 'overflow-x:scroll;',
                          DT::dataTableOutput("s_summary_table")
                      )
             ),
             tabPanel("PCA",
                      div(style = 'overflow-x:scroll;',
                          plotOutput("s_pcabiplot"),
                          plotOutput("s_pcascreeplot")
                      )
             ),
             tabPanel("QC Box Plots",
                      div(style = 'overflow-x:scroll;',
                          plotlyOutput("s_global_boxplot")
                      )
             )
      )
    } else {
      values$selected_dataset_name<-"No dataset selected"
      div("View and explore a simulated dataset by clicking on the tile on the heatmap that corresponds to the requested number of proteins and biological replicates. If required, vary the range of datasets to be simulated by modifying the parameters.")
    }
    
  })
  
  ## L2 Explore simulated data / Specified dataset / Data tables
  output$s_prot_abundance_table = DT::renderDataTable({
    values$s_prot_abundance
  })
  output$s_sample_annotation_table = DT::renderDataTable({
    values$s_groups
  })
  output$s_summary_table = DT::renderDataTable({
    values$s_summary
  })
  
  ## L2 Explore simulated data / Specified dataset / Plot: simulated_grid
  output$simulated_grid<-renderPlotly({
    protein_desc<-1/input$n_protein
    nproteins <- length(unique(values$quant_data$ProcessedData$PROTEIN))
    m_prot <- floor(nproteins*protein_desc)
    ngroup <- length(unique(values$quant_data$ProcessedData$GROUP_ORIGINAL))
    sample_incr <- input$sample_incr
    n_sample <- input$n_sample
    train_size <- seq.int(from = sample_incr, to = sample_incr * n_sample, length.out = n_sample)
    train_size <- train_size * ngroup
    protein_num <- seq.int(from = m_prot, to = nproteins, by = m_prot)
    
    rownames <- as.character(train_size)
    colnames <- as.character(protein_num)
    
    simulated_grid<-matrix(NA, nrow=length(rownames), ncol=length(colnames), dimnames = list(rownames,colnames))
    for (m in rownames(simulated_grid)) {
      for (n in colnames(simulated_grid)) {
        simulated_grid[m,n] <- as.numeric(m)*as.numeric(n) 
      }
    }
    plot_ly(
      x=colnames(simulated_grid),
      y=rownames(simulated_grid),
      z=simulated_grid,
      type="heatmap",
      source="simulated_grid_plot"
    ) %>%
      layout(xaxis=list(title="Number of Proteins", type="category"), yaxis=list(title="Number of Biological Replicates", type="category"))
  })
  
  # T4 Experiment Simulation
  ## Main Content
  output$analyse_simulated_content<-renderUI({
    if (is.null(values$is_imported)) {
      tags$b("Please use the 'Import Data' menu to import a proteome dataset.")
    } else {
      tagList(
        fluidRow(
          box(width = 9, title = values$simulated_classification_title,
              htmlOutput("simulated_classification_heatmap")
              ),
          box(title="Train Models", width = 3,
              radioButtons("train_type", label="Select type of analysis", 
                           choices=list("Between datasets (Single classifier)"="between_datasets", 
                                        "Between classifiers (Single dataset)"="between_classifiers")
              ),
              htmlOutput("train_options"),
              numericInput("iter", label="Iteration Count", value=10),
              actionButton("start", label="Begin Simulation")
          )
        ),
        fluidRow(
          htmlOutput("inspect_model")
        )
      )
    }
  })
  
  ## L2 Experiment Simulation / Training Options
  output$train_options <- renderUI({
    if (!is.null(input$train_type)) {
      switch(input$train_type,
             "between_datasets" = {
               # numericInput("n_sample", "Number of proteins", value=input$s_n_sample),
               # 
               # n_sample = input$n_sample,
               # sample_incr = input$sample_incr,
               # protein_desc = 1/input$protein_desc,
               selectInput("classifier", "Select Classifier", choices=list("Random Forest" = "rf", "SVM" = "svmLinear", "Naive Bayes" = "naive_bayes", "Partial Least Squares" = "pls", "Neural Net (3 layers)" = "nnet"))
             },
             "between_classifiers" = {
               values$chosen_dataset <- event_data("plotly_click", source = "model_plot")
               if (length(values$chosen_dataset)) {
                 values$between_classifiers_prot <- values$chosen_dataset[["x"]]
                 values$between_classifiers_samples <- values$chosen_dataset[["y"]]
               }
               tagList(
                 tags$b("Choose a dataset by clicking on its tile on the heatmap"),
                 actionButton("between_classifiers_declare", "Select Dataset"),
                 tags$p(paste("Selected dataset:",values$between_classifiers_prot," Proteins, ",values$between_classifiers_samples," BioReplicates"))
               )
             }
             )
    }
  })

  ## L2 Experiment Simulation / Heatmap
  output$simulated_classification_heatmap<-renderUI({
    if (is.null(values$simulated_classification_status)) {
      values$simulated_classification_title = "Waiting for trained models..."
      tagList(
        div('Train classifiers on generated fake datasets by clicking "Begin Simulation" on the right. Optionally, change the classifier to be used, or choose to analyse cross-classifier performance for a single dataset on the right.'),
        textOutput("main_output")
      )
    } else {
      values$simulated_classification_title = "Trained Models"
      tabBox(id="simulated_classification_wrapper", width="100%",
             tabPanel("Heatmap",
                      plotlyOutput("model_heatmap")
             ),
             tabPanel("Line Graphs",
                      htmlOutput("lineplot_xvar_selector"),
                      htmlOutput("lineplot_trendline_selector"),
                      plotOutput("lineplot") 
             )
      )
    }
  })

  # L2 Experiment Simulation / Lineplot Type Selectors
  output$lineplot_xvar_selector <- renderUI({
    selectInput("lineplot_xvar", "X variable:", choices = c("Protein number", "Sample size"))
  })
  output$lineplot_trendline_selector <- renderUI({
    selectInput("lineplot_trendline", "Trendline Type", choices = c("None", "Linear", "Polynomial"))
  })
  
  # L1 Inspect Model
  output$inspect_model<-renderUI({
    values$chosen_dataset <- event_data("plotly_click", source="model_plot")
    if (length(values$chosen_dataset)) {
      
      chosen_prot <- as.numeric(values$chosen_dataset[["x"]])
      chosen_rep <- as.numeric(values$chosen_dataset[["y"]])
      message("2 Done")
      
      prot_index <- which(values$result$paramlog$protein_num==chosen_prot)
      
      train_index <- which(values$result$paramlog$train_size==chosen_rep)
      message("3 Done")
      message(prot_index, train_index)
      # Get model with the highest accuracy from all iterations on the given dataset
      
      iter <- length(values$result$results)
      accuracy <- vector()
      for (i in 1:iter) {
        message(i, " Model Checked")
        currentmodel <- values$result$results[[i]][[prot_index]][[train_index]][[2]]
        currentmodel.pred <- predict(currentmodel, values$valid_x) #Predict validation data
        accuracy[i] <- sum(diag(table(currentmodel.pred, values$valid_y))) / length(currentmodel.pred)
        highestacc <- which(accuracy==max(accuracy))[[1]]
      }
      message("highestacc: ", highestacc)
      message("4 Done")
      
      values$chosen_model <- values$result$results[[highestacc]][[prot_index]][[train_index]][[2]]
      values$chosen_model.pred <- predict(values$chosen_model, values$valid_x)
      values$confusionMatrix <- confusionMatrix(data = as.factor(as.numeric(values$chosen_model.pred)), reference = as.factor(values$valid_y))
      message("5 Done")
      
      d1 <- as.data.frame(values$confusionMatrix$overall)
      colnames(d1)<-"Value"
      d2 <- as.data.frame(values$confusionMatrix$byClass)
      colnames(d2)<-"Value"
      # c2 <- rbind(d1, d2)
      # c1 <- rownames(c2)
      values$chosen_model.metrics <- rbind(d1, d2)
      message("6 Done")
      
      tabBox(id="inspect_model_wrapper", width=12,
             tabPanel("Confusion Matrix",
                      div(style="display: inline-block;vertical-align:top; width: 75%;",
                          plotOutput("confusion_matrix")
                          ),
                      div(style="display: inline-block;vertical-align:top; width: 24%;",
                          DT::dataTableOutput("metrics_table")                      
                          )
                      ),
             tabPanel("Predictive Features",
                      plotOutput("predictive_features")
                      ),
             tabPanel("LIME",
                      pickerInput("test_set_picker", 
                                  label="Select observations to test model with", 
                                  choices = as.character(rownames(values$valid_x)),
                                  selected=1,
                                  options = list(
                                    `actions-box` = TRUE, 
                                    size = 10
                                    ),
                                  multiple = TRUE),
                      sliderInput("feature_picker",
                                   label="Number of features to display",
                                   value=4,
                                   min=1,
                                   max=10,
                                   step=1,
                                   round=TRUE),
                      plotOutput("lime_plot", height="auto")
                      )
             )
    }
  })
  
  # L2 Inspect model / LIME
  output$lime_plot <- renderPlot({
    test_set <- as.integer(input$test_set_picker)
    explainer <- lime(values$valid_x[-test_set,], as_classifier(values$chosen_model), bin_continuous = TRUE, quantile_bins = FALSE)
    explanation <- explain(values$valid_x[test_set, ], explainer, n_labels = 1, n_features = input$feature_picker)
    plot_features(explanation, ncol = 1)
  }, height = function(){
    (70 + input$feature_picker * 50) * length(input$test_set_picker)
  })
  
  
  # Plots & Tables / confusion_matrix
  output$confusion_matrix <- renderPlot(
    ggplot(as.data.frame(values$confusionMatrix$table), aes(x=Prediction, y=Reference, fill=Freq)) +
      geom_tile() +
      geom_text(aes(label=Freq))
  )
  
  # Plots & Tables / metrics 
  output$metrics_table = DT::renderDataTable({
    values$chosen_model.metrics
  })
  
  # Plots & Tables / variable importance
  output$predictive_features <- renderPlot({
    plot(varImp(values$chosen_model), top=10)
  })
  

  # Actual Simulation Function
  observeEvent(input$start, {
    withCallingHandlers({
      shinyjs::html("text", "")
      shinyjs::html(id = "main_output", html = "")
      # shinyjs::runjs(
      #   var element = document.getElementById("mainOutput");
      #   element.scrollTo = element.scrollHeight;
      # )
      values$result<-ParSimulatedClassificationCrossDatasets(values$quant_data,
                                             values$parameters,
                                             n_sample = input$n_sample,
                                             sample_incr = input$sample_incr,
                                             n_protein = input$n_protein,
                                             iter = input$iter,
                                             classifier = input$classifier,
                                             use_caret = TRUE)
      rownames(values$result$meanPA) <- gsub("prot", "", rownames(values$result$meanPA))
      colnames(values$result$meanPA) <- gsub("tra", "", colnames(values$result$meanPA))
      values$gg_meanPA<-melt(values$result$meanPA)
      
      # Clean up the validation set as output from parameters
      values$valid_x <- as.data.frame(values$parameters$X)
      rownames(values$valid_x)<-as.numeric(as.factor(rownames(values$valid_x)))              
      values$valid_y <- as.factor(values$parameters$Y)
      values$valid_y <- as.numeric(values$valid_y)
      message("1 Done")
      
      values$simulated_classification_status <- "completed"
      
    },
    message = function(m) {
      shinyjs::html(id = "main_output", html = paste(m$message,"<br>", sep=" "), add = TRUE)
    })
  })
  
  ### For Simulated
  ## PCA
  output$s_pcabiplot <- renderPlot({
    ggbiplot(values$s_pca, ellipse=TRUE, var.axes=FALSE, groups=values$s_sample_annotation)
    # plot_ly(values$pca, x=Comp.1, y=Comp.2, text=rownames(values$pca), mode="markers", marker=list(size=11))
  })
  
  output$s_pcascreeplot <- renderPlot({
    ggscreeplot(values$s_pca)
  })
  
  ## Global QC Box Plot
  output$s_global_boxplot <- renderPlotly({
    plot_ly(values$s_quant_data, y=~LogIntensities, x=~as.character(originalRUN), color=~Group, type = "box") %>%
      layout(xaxis=list(title="Biological Replicate"), yaxis=list(title="Log Intensities"))
    
  })
  
  # Experiment Plots
  output$model_heatmap<-renderPlotly({
    plot_ly(
      y=colnames(values$result$meanPA),
      x=rownames(values$result$meanPA),
      z=values$result$meanPA,
      type="heatmap",
      source="model_plot"
    ) %>%
      layout(xaxis=list(title="Number of Proteins", type="category"), yaxis=list(title="Number of Biological Replicates", type="category"))
    # 
    # hm.palette <- colorRampPalette(rev(brewer.pal(9, 'RdBu')), space='Lab')
    # ggplot(data=values$gg_meanPA, aes(x=Var1, y=Var2, fill=value)) +
    #  geom_tile() +
    #  # geom_text(aes(label=round(value,3))) +
    #  coord_fixed(ratio=max(values$gg_meanPA$Var1, na.rm = TRUE)/max(values$gg_meanPA$Var2, na.rm = TRUE)) +
    #  scale_fill_viridis() +
    #  labs(title="Heatmap", x="Protein number", y='Mean accuracy')

    # plot_ly(z=~values$result$meanPA, type="heatmap")
  })
  
  output$lineplot<-renderPlot({
      ## ggplot needs a long format dataframe
      ## get the mean accuracy
      meandata <- melt(values$result$meanPA)
      colnames(meandata) <- c("protein_num", "sample_size", "mean_acc")
      
      switch(input$lineplot_xvar,
        "Protein number" =
          p <- ggplot(data = meandata, aes(x= protein_num, y = mean_acc, group = sample_size, colour = sample_size)) +
          geom_point() +
          # geom_line() +
          # geom_smooth(method = input$smoothing_method) +
          labs(title="Sample size estimation", x="Protein number", y="Mean accuracy") +
          guides(color=guide_legend(title="Sample size")),
        "Sample size" =
          p <- ggplot(data = meandata, aes(x= sample_size, y = mean_acc, group = protein_num, colour = protein_num)) +
          geom_point() +
          # geom_smooth(method = input$smoothing_method) +
          labs(title="Sample size estimation", x="Sample size", y='Mean accuracy') +
          guides(color=guide_legend(title="Protein number")),
        )
      
      switch(input$lineplot_trendline,
          "None" =
            p <- p + geom_line(),
          "Linear" =
            p <- p + geom_smooth(method = lm, se = FALSE),
          "Polynomial" =
            p <- p + geom_smooth(method = lm, formula = y ~ poly(x, 3), se = FALSE)
      )

      ## make the plot
      p + geom_line()
      return(p)

  })
}