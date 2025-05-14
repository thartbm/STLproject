library(Reach)
library(optimx)
library(dplyr)
library(ggplot2)
library(tidyr)

Download_data <- function() {
  Reach::downloadOSFdata(repository='6m24e',
                         filelist=list('data'=c('data.zip')), 
                         folder='data', 
                         unzip=TRUE, 
                         removezips=TRUE)
}

# Preprocess participant data
Preprocess_file <- function(id, condition) {
  # Read file for participant id
  filename <- sprintf('data/%s_performance.csv', id) 
  df <- read.csv(filename, stringsAsFactors = FALSE)
  
  # Extract relevant data 
  idx <- which(df$label == condition)
  response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
  rotation <- df$rotation[idx]
  
  # Exclude outliers and filter corresponding elements
  filtered_idx <- idx[response >= -60 & response <= 60]
  filtered_response <- response[response >= -60 & response <= 60]
  filtered_rotation <- rotation[response >= -60 & response <= 60]
  
  # Store into data frame
  output <- data.frame(rotation = filtered_rotation, 
                       response = filtered_response, 
                       trial = filtered_idx)
  output$participant <- id
  
  # Normalize data
  output$response[which(output$rotation > 0)] <- -1 * output$response[which(output$rotation > 0)]
  output$rotation[which(output$rotation < 0)] <- -1 * output$rotation[which(output$rotation < 0)]
  
  return(output)
}

# Function to preprocess data and write to CSV
RawData_CSV <- function(data_portion = "all", file_name) {
  # Initialize a data frame to store all processed data
  all_data <- data.frame()
  
  # Setup a vector of id strings
  csv_files <- list.files("data", pattern = "_performance.csv", full.names = TRUE)
  
  # Loop through each CSV file and preprocess data
  for (file in csv_files) {
    # Extract the identifier (xxxxxx) from the file name
    identifier <- gsub("_performance.csv", "", basename(file))
    
    # Read file for participant id
    df <- read.csv(file, stringsAsFactors = FALSE)
    
    # Extract relevant data
    idx <- which(df$label == 'stl-arc-rotation')
    response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
    rotation <- df$rotation[idx]
    
    # Exclude outliers and filter corresponding elements
    filtered_idx <- idx[response >= -60 & response <= 60]
    filtered_response <- response[response >= -60 & response <= 60]
    filtered_rotation <- rotation[response >= -60 & response <= 60]
    
    # Store into data frame
    output <- data.frame(Participant_ID = rep(identifier, length(filtered_idx)),
                         Rotation = filtered_rotation,
                         Response = filtered_response)
    
    # Normalize data
    output$Response[output$Rotation > 0] <- -1 * output$Response[output$Rotation > 0]
    output$Rotation[output$Rotation < 0] <- -1 * output$Rotation[output$Rotation < 0]
    
    # Debug print to check the number of rows before subsetting
    cat("Total rows for", identifier, "before subsetting:", nrow(output), "\n")
    
    # Subset the data based on the specified portion
    if (data_portion == "first_half") {
      output <- output[1:floor(nrow(output) / 2), ]
    } else if (data_portion == "last_half") {
      output <- output[(floor(nrow(output) / 2) + 1):nrow(output), ]
    }
    
    # Debug print to check the number of rows after subsetting
    cat("Total rows for", identifier, "after subsetting:", nrow(output), "\n")
    
    # Add data for the current participant to the all_data data frame
    all_data <- rbind(all_data, output)
  }
  
  # Write processed data to CSV file
  write.csv(all_data, file_name, row.names = FALSE)
}


getidentifiers <- function() {
  
  # Setup a vector of id strings
  csv_files <- list.files("data", pattern = "_performance.csv", full.names = TRUE)
  
  # Create an empty list to store identifiers
  identifiers_list <- c()
  
  # Loop through each CSV file and extract the identifier
  for (file in csv_files) {
    # Extract the identifier (xxxxxx) from the file name
    identifier <- gsub("_performance.csv", "", basename(file))
    
    # Add the identifier to the list
    identifiers_list <- c(identifiers_list, identifier)
  }
  
  return(identifiers_list)
  
}


Preprocess_all_files <- function() {
  
  triallabels <- c('stl-target-rotation', 'stl-arc-rotation')
  identifiers <- getidentifiers()
  for (condition in triallabels) {
    
    condition_data <- NA
    
    for (identifier in identifiers) {
      processed_data <- Preprocess_file(identifier, condition)
      if (is.data.frame(condition_data)) {
        condition_data <- rbind(condition_data, processed_data)
      } else {
        condition_data <- processed_data
      }
    }
    
    write.csv( x= condition_data,
               file= sprintf('data/%s.csv', condition),
               quote=TRUE,
               row.names =FALSE)
    
    avg_condition_data <- aggregate(response ~ rotation + participant, data=condition_data, FUN=mean, rm.na=TRUE)
    write.csv( x         = avg_condition_data,
               file      = sprintf('data/avg_%s.csv', condition),
               quote     = TRUE,
               row.names = FALSE)
    
  }
  
}

# STL Prediction Function
STLPredict <- function(par, rotations) {
  
  return((rotations * par['s']) * (dnorm(rotations, mean=0, sd= par['w']) / dnorm(0, mean=0, sd= par['w'])))
}

# Function to establish & return predicted adaptions for given fit data and rotations  
STLPredictions <- function(fits_data, rotations) {
  results <- data.frame(participant = character(), rotations = numeric(), predictions = numeric(), stringsAsFactors = FALSE)
  
  fits_data$slope <- as.numeric(fits_data$slope)
  fits_data$width <- as.numeric(fits_data$width)
  
  for (i in 1:nrow(fits_data)) {
    participant <- fits_data$participant[i]
    par <- c(s = fits_data$slope[i], w = fits_data$width[i])
    predictions <- STLPredict(par = par, rotations = rotations)
    
    participant_results <- data.frame(participant = rep(participant, length(rotations)),
                                      rotations = rotations,
                                      predictions = predictions)
    
    results <- rbind(results, participant_results)
  }
  
  return(results)
}    


# Function to establish & return predicted adaptions for given fit data and rotations  
STLConditionPredictions <- function(fits_data, rotations, condition) {
  # Filter fits_data based on the condition
  fits_data <- fits_data[fits_data$condition == condition, ]
  
  # Convert slope and width columns to numeric
  fits_data$slope <- as.numeric(fits_data$slope)
  fits_data$width <- as.numeric(fits_data$width)
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each row of the fits_data
  for (i in seq_len(nrow(fits_data))) {
    participant <- fits_data$participant[i]
    par <- c(s = fits_data$slope[i], w = fits_data$width[i])
    
    # Get predictions for the current participant
    predictions <- STLPredict(par = par, rotations = rotations)
    
    # Store results in a temporary data frame
    participant_results <- data.frame(
      participant = rep(participant, length(rotations)),
      rotations = rotations,
      predictions = predictions
    )
    
    # Append participant results to the results list
    results_list[[i]] <- participant_results
  }
  
  # Combine all participant results into a single data frame
  results <- do.call(rbind, results_list)
  
  return(results)
}

# STL Error Function
STLErrors <- function(par, rotations, deviations) {
  
  # Generate model predictions using STLpredict function
  model_predictions <- STLPredict(par, rotations)
  
  # Calculate squared errors between model predictions and actual deviations
  squared_errors <- (deviations - model_predictions)^2
  
  # Compute the mean squared error (MSE)
  MSE <- mean(squared_errors)
  
  return(MSE)
}

# STL Grid Search Function
STLGridsearch <- function(rotations, deviations) {
  # Define parameter grids
  s_values <- seq(0.1, 1, length.out = 10)  # Example range for parameter s
  w_values <- seq(1, 60, length.out = 60)  # Example range for parameter w
  
  # Generate all combinations of parameters
  parameter_combinations <- expand.grid(s = s_values, w = w_values)
  
  MSE <- apply(parameter_combinations, FUN=STLErrors, MARGIN=c(1), rotations = rotations, deviations = deviations)
  
  top_parameters <- parameter_combinations[order(MSE)[1:10], ]
  
  return(top_parameters)
}


STLOptimization <- function(top_parameters, rotations, deviations) { # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( top_parameters,
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=STLErrors,
                            method = "L-BFGS-B", 
                            lower = c(0.1, 1), 
                            upper = c(1, 60),
                            rotations = rotations,
                            deviations = deviations
                     ))
  
  win <- allfits[order(allfits$value)[1], ]
  return(unlist(win[1:2]))
}

# Returns individual STL Fit parameters for a list of participants
STLIndividualFits <- function(data) {
  # Create an empty data frame to store the best fits for all participants
  all_fits_df <- data.frame(participant = character(),
                            slope = numeric(),
                            width = numeric(),
                            stringsAsFactors = FALSE)
  
  identifiers <- unique(data$participant)
  
  # Loop through each participant
  for (identifier in identifiers) {
    
    cat(sprintf('fitting: %s\n', identifier))
    
    # get row indices for participant
    idx <- which(data$participant == identifier)
    
    # Extract rotations and deviations
    rotations <- data$rotation[idx]
    deviations <- data$response[idx]
    
    # Run grid search to find best parameters
    top_parameters <- STLGridsearch(rotations, deviations)
    
    # Run optimization on the best starting positions
    winning_par <- STLOptimization(top_parameters, rotations, deviations)
    
    # Add the participant id, slope, and width to the data frame
    all_fits_df[nrow(all_fits_df) + 1, ] <- c(identifier, winning_par[1], winning_par[2])
  }
  
  # Set column names
  colnames(all_fits_df) <- c("participant", "slope", "width")
  
  return(all_fits_df)
  
}

STLindividualConditionFits <- function() {
  
  conditions <- c('target', 'arc')
  
  all_fits <- NA
  
  for (condition in conditions) {
    
    data <- read.csv(sprintf('data/stl-%s-rotation.csv', condition), stringsAsFactors = FALSE)
    
    fits <- STLIndividualFits(data = data)
    fits$condition <- condition
    
    if (is.data.frame(all_fits)) {
      all_fits <- rbind(all_fits, fits)
    } else {
      all_fits <- fits
    }
    
  }
  
  write.csv(all_fits, 'data/condition_fits.csv', quote=FALSE, row.names=FALSE)
  
}


# Returns a single (or weighted) set of STL Fit parameters for a set of participant data
STLFit <- function(participant, data) {
  # Read raw data from CSV
  #data <- read.csv(data)
  
  # If participant is a single ID, convert it to a list
  if (!is.vector(participant)) {
    participant <- vector(participant)
  }
  
  # Extract data for all participants
  all_data <- NA
  for (id in participant) {
    if (is.data.frame(all_data)) {
      all_data <- rbind(all_data, data[which(data$participant == id),])
    } else {
      all_data <- data[which(data$participant == id),]
    }
  }
  
  # Extract rotations and deviations from combined data
  rotations <- all_data$rotation
  deviations <- all_data$response
  
  # Run grid search to find best parameters
  top_parameters <- STLGridsearch(rotations, deviations)
  
  # Run optimization on the best parameters
  winning_par <- STLOptimization(top_parameters, rotations, deviations)
  
  # Create a data frame for the results
  results_df <- data.frame(Participant_ID = "Combined",
                           Slope = winning_par[1],
                           Width = winning_par[2],
                           stringsAsFactors = FALSE)
  
  return(results_df)
}

# Baselines and derives Exponential Model Fits from participant Fixed Rotation data 
STLExponentialFits <- function(ids) {
  if (!is.vector(ids)) {
    ids <- c(ids)  # Ensure ids is a vector
  }
  
  results <- data.frame(participant = character(), lambda = numeric(), N0 = numeric(), stringsAsFactors = FALSE)
  
  for (id in ids) {
    filename <- sprintf('data/%s_performance.csv', id)
    df <- read.csv(filename, stringsAsFactors = FALSE)
    
    # Extract baseline data
    bias <- median(df$reachdeviation_deg[which(df$label == 'fixed-rotation-baseline')[c(11:60)]])
    learning <- df$reachdeviation_deg[which(df$label == 'fixed-rotation')] - bias
    
    # Derive Exponential Fits
    fit <- Reach::exponentialFit(signal = learning)
    
    # Initialize a data frame to store Exponential Fit data
    participant_result <- data.frame(participant = id, lambda = fit['lambda'], N0 = fit['N0'])
    
    # Store Exponential Fit Data
    results <- rbind(results, participant_result)
  }
  
  return(results)
}


# Obtains graph-able Exponential Model values for participants
STLExponentialModel <- function(participants_results, mode = 'learning', setN0 = NULL, points) {
  
  participants_results$lambda <- as.numeric(participants_results$lambda)
  participants_results$N0 <- as.numeric(participants_results$N0)
  
  # Initialize an empty data frame to hold all the plot data
  all_plot_data <- data.frame()
  
  # Iterate over each participant's results
  for (i in 1:nrow(participants_results)) {
    participant_result <- participants_results[i, ]
    lambda <- participant_result$lambda
    N0 <- participant_result$N0
    participant <- participant_result$participant
    
    # Generate the time points (for example, 0 to 100)
    time_points = seq(0, points, by = 1)
    
    # Generate the model values using the exponentialModel function
    model_output <- Reach::exponentialModel(par = c('lambda' = lambda, 'N0' = N0), timepoints = time_points, mode = mode, setN0 = setN0)
    
    # Create a data frame for the current participant's plot data
    plot_data <- data.frame(time = model_output$trial, value = model_output$output, participant = participant)
    
    # Combine with the overall plot data
    all_plot_data <- rbind(all_plot_data, plot_data)
  }
  return(all_plot_data)
}

PlotReachdata <- function(identifiers_list, filename, predictions_data, condition) {
  
  dir_path <- "C:/Users/prahi/Desktop/ErrorsSTL-MSc/ErrorsSTL/doc"
  
  pdf_name <- paste0(dir_path, "/", filename, ".pdf")
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  pdf(pdf_name, width = 8.5, height = 11)  # Open PDF device
  
  # Initialize a list to store plots
  plot_list <- list()
  
  # Loop through each participant
  for (identifier in identifiers_list) {
    # Extract data for the current participant using the Preprocess_file function
    participant_data <- Preprocess_file(identifier, condition)
    
    # Calculate Mean Deviation per rotation value
    mean_deviation <- participant_data %>%
      group_by(rotation) %>%
      summarize(mean_deviation = mean(response))
    
    # Extract predicted data for the current participant
    participant_predictions <- predictions_data %>%
      filter(participant == identifier)
    
    p <- ggplot(participant_data, aes(x = rotation, y = response)) +
      geom_point(color = "black") +
      
      # Plot mean data points with red line
      geom_line(data = mean_deviation, aes(x = rotation, y = mean_deviation), color = "red", alpha = 0.6) +
      
      # Plot predicted data points with blue line
      geom_line(data = participant_predictions, aes(x = rotations, y = participant_predictions[[3]]), color = "blue", alpha = 0.6) +
      
      # Add a dotted grey line at y = 0
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +    
      
      labs(title = identifier,
           x = "Rotation in Degrees(°)", 
           y = "Deviation in Degrees(°)") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    # Add the plot to the list
    plot_list <- c(plot_list, list(p))
    
    # If we have 6 plots, print them to a single page
    if (length(plot_list) == 6) {
      do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 2))
      plot_list <- list()  # Reset the list
    }
  }
  
  # If there are remaining plots not yet printed, print them on the last page
  if (length(plot_list) > 0) {
    do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 2))
  }
  
  dev.off()  # Close PDF device
}

# Function to plot exponential model values 
STLExponentialPlot <- function(data, rot, vert_line = NULL, trials = NULL){
  
  main_title <- paste(paste0(rot,"°"), "Adaptation Over", trials, "Trials")
  
  # Allows to plot for a certain number of trials
  if (is.null(trials)){
    
    ggplot(data, aes(x = time, y = value, color = participant)) +
      geom_line(alpha = 0.7) +
      labs(title = main_title,
           x = "Trials",
           y = "Model Value") +
      theme_minimal() +
      guides(color = "none") +
      geom_vline(xintercept = vert_line, linetype = "dotted", color = "black")
    
    
  } else{
    # Add the trials column, resetting at 0 for each new participant
    data <- data %>%
      group_by(participant) %>%
      mutate(trials = row_number() - 1) %>%
      ungroup()
    
    # Filter data to include only up to the specified number of trials
    data <- data %>%
      filter(trials <= trials)
    
    ggplot(data, aes(x = trials, y = value, color = participant)) +
      geom_line(alpha = 0.7) +
      labs(title = main_title,
           x = "Trials",
           y = "Model Value") +
      theme_minimal() +
      guides(color = "none") +
      geom_vline(xintercept = vert_line, linetype = "dotted", color = "black") +
      scale_x_continuous(limits = c(0, trials))
    
  }
}

# Function to plot predicted STL adaptation values
STLPredictPlot <- function(data, vert_line = NULL, main){
  
  p <- ggplot(data, aes(x = rotations, y = predictions, group = participant)) +
    geom_line(color = "red", alpha = 0.6) +  # Lower alpha for more transparency
    labs(title = main,
         x = "Rotations in Degrees(°)",
         y = "Predicted Adaptation in Degrees(°)") +
    theme_minimal()
  
  # Add vertical line if specified
  if (!is.null(vert_line)) {
    p <- p + geom_vline(xintercept = vert_line, linetype = "dotted", color = "black")
  }
  
  print(p)  # Ensure the plot is printed
}

  # Function to plot Arc-Circle Plots
PlotConditionPredictions <- function(identifiers_list, filename, predictions_data) {
  
  dir_path <- "C:/Users/prahi/Desktop/ErrorsSTL-MSc/ErrorsSTL/doc/"
  
  pdf_name <- paste0(dir_path, filename, "_Reach_Data.pdf")
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  pdf(pdf_name, width = 8.5, height = 11)  # Open PDF device
  
  
  data_file <- read.csv("C:/Users/prahi/Desktop/ErrorsSTL-MSc/ErrorsSTL/raw_arc.csv", stringsAsFactors = FALSE)
  
  
  # Initialize a list to store plots
  plot_list <- list()

  # Loop through each participant
  for (identifier in identifiers_list) {
    # Extract data for the current participant from the data file
    participant_data <- data_file[data_file$Participant_ID == identifier, ]
  
    # Extract rotations and deviations
    rotations <- participant_data$Rotation
    deviations <- participant_data$Response
  
    # Calculate Mean Deviation per rotation value
    mean_deviation <- participant_data %>%
      group_by(Rotation) %>%
      summarize(mean_deviation = mean(Response))
  
    # Extract predicted data for the current participant
    participant_predictions <- predictions_data %>%
      filter(participant == identifier)
  
    p <- ggplot(participant_data, aes(x = Rotation, y = Response)) +
    
      # Plot predicted data points with blue line
      geom_line(data = participant_predictions, aes(x = rotations, y = predictions_arc), color = "blue", alpha = 0.6) +
      # Plot predicted data points with blue line
      geom_line(data = participant_predictions, aes(x = rotations, y = predictions_circle), color = "red", alpha = 0.6) +
      # Plot the difference with green line
      geom_line(data = participant_predictions, aes(x = rotations, y = predictions_diff), color = "green", alpha = 0.6) +
    
      # Add a dotted grey line at y = 0
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey") + 
    
      labs(title = identifier,
           x = "Rotation in Degrees(°)",
           y = "Deviation in Degrees(°)") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
    # Add the plot to the list
    plot_list <- c(plot_list, list(p))
  
    # If we have 6 plots, print them to a single page
    if (length(plot_list) == 6) {
      do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 2))
      plot_list <- list()  # Reset the list
    }
  }

  # If there are remaining plots not yet printed, print them on the last page
  if (length(plot_list) > 0) {
    do.call(grid.arrange, c(plot_list, ncol = 3, nrow = 2))
  }

  dev.off()  # Close PDF device
  
}

STLscatter <- function(data, xax, yax, col, by, title){
  dir_path <- "C:/Users/prahi/Desktop/ErrorsSTL-MSc/ErrorsSTL/doc"
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  pdf_name <- paste0(dir_path, "/", title, ".pdf")
  
  # Reshape the data based on the `by` argument
  # Reshape the data based on the `by` argument
  reshaped_data <- data %>%
    select(participant, !!sym(by), condition) %>%
    pivot_wider(names_from = condition, values_from = !!sym(by))
  
  # Define the plot
  plot <- ggplot(reshaped_data, aes(x = !!sym(xax), y = !!sym(yax), color = "black")) +
    geom_point(color ="black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey") +
    labs(title = title,
         x = xax,
         y = yax) +
    theme_minimal()
  
  # Save the plot as a PDF
  ggsave(pdf_name, plot = plot, width = 8.5, height = 11)
}

STLPlotPredictions <- function(data, vert_line = NULL, xax, yax, col, title) {
  dir_path <- "C:/Users/prahi/Desktop/ErrorsSTL-MSc/ErrorsSTL/doc"
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  pdf_name <- paste0(dir_path, "/", title, ".pdf")
  
  # Define the plot
  plot <- ggplot(data, aes(x = !!sym(xax), y = !!sym(yax), group = participant)) +
    geom_line(color = col) +
    labs(title = title,
         x = xax,
         y = yax) +
    theme_minimal() +
    geom_vline(xintercept = vert_line, linetype = "dotted", color = "black")  
  
  # Save the plot as a PDF
  ggsave(pdf_name, plot = plot, width = 8.5, height = 11)
}


