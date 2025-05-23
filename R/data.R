library(Reach)
library(optimx)

Download_data <- function() {
  Reach::downloadOSFdata(repository='6m24e',
                         filelist=list('data'=c('data.zip', 'demographics.csv')), 
                         folder='data', 
                         unzip=TRUE, 
                         removezips=TRUE)
}

# pre-processing data -----

# Preprocess participant data
Preprocess_file <- function(id, condition) {
  # Read file for participant id
  filename <- sprintf('data/%s_performance.csv', id) 
  df <- read.csv(filename, stringsAsFactors = FALSE)
  
  # add super block numbers?
  df$superblock <- 1
  labidx <- which(df$label == 'stl-baseline')
  diffidx <- diff(labidx) > 1
  
  incidx <- which( df$label == 'stl-baseline' )[which(diff(which(df$label == 'stl-baseline')) > 1) + 1]
  for (ii in incidx) {
    df$superblock[c(ii:dim(df)[1])] <- df$superblock[c(ii:dim(df)[1])] + 1
  }
  
  # Extract relevant data
  idx <- which(df$label == condition)
  response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
  rotation <- df$rotation[idx]
  superblock <- df$superblock[idx]
  
  # Exclude outliers and filter corresponding elements
  filtered_idx <- idx[response >= -60 & response <= 60]
  filtered_response <- response[response >= -60 & response <= 60]
  filtered_rotation <- rotation[response >= -60 & response <= 60]
  superblock <- superblock[response >= -60 & response <= 60]
  
  # Store into data frame
  output <- data.frame(rotation = filtered_rotation, 
                       response = filtered_response,
                       superblock = superblock)
  output$participant <- id
  
  # Normalize data
  output$response[which(output$rotation > 0)] <- -1 * output$response[which(output$rotation > 0)]
  output$rotation[which(output$rotation < 0)] <- -1 * output$rotation[which(output$rotation < 0)]
  
  return(output)
}

# # Function to preprocess data and write to CSV
# RawData_CSV <- function(data_portion = "all", file_name) {
#   # Initialize a data frame to store all processed data
#   all_data <- data.frame()
#   
#   # Setup a vector of id strings
#   csv_files <- list.files("data", pattern = "_performance.csv", full.names = TRUE)
#   
#   # Loop through each CSV file and preprocess data
#   for (file in csv_files) {
#     # Extract the identifier (xxxxxx) from the file name
#     identifier <- gsub("_performance.csv", "", basename(file))
#     
#     # Read file for participant id
#     df <- read.csv(file, stringsAsFactors = FALSE)
#     
#     # Extract relevant data
#     idx <- which(df$label == 'stl-target-rotation')
#     response <- df$reachdeviation_deg[idx + 1] - df$reachdeviation_deg[idx - 1]
#     rotation <- df$rotation[idx]
#     
#     # Exclude outliers and filter corresponding elements
#     filtered_idx <- idx[response >= -60 & response <= 60]
#     filtered_response <- response[response >= -60 & response <= 60]
#     filtered_rotation <- rotation[response >= -60 & response <= 60]
#     
#     # Store into data frame
#     output <- data.frame(Participant_ID = rep(identifier, length(filtered_idx)),
#                          Rotation = filtered_rotation,
#                          Response = filtered_response)
#     
#     # Normalize data
#     output$Response[output$Rotation > 0] <- -1 * output$Response[output$Rotation > 0]
#     output$Rotation[output$Rotation < 0] <- -1 * output$Rotation[output$Rotation < 0]
#     
#     # Debug print to check the number of rows before subsetting
#     cat("Total rows for", identifier, "before subsetting:", nrow(output), "\n")
#     
#     # Subset the data based on the specified portion
#     if (data_portion == "first_half") {
#       output <- output[1:floor(nrow(output) / 2), ]
#     } else if (data_portion == "last_half") {
#       output <- output[(floor(nrow(output) / 2) + 1):nrow(output), ]
#     }
#     
#     # Debug print to check the number of rows after subsetting
#     cat("Total rows for", identifier, "after subsetting:", nrow(output), "\n")
#     
#     # Add data for the current participant to the all_data data frame
#     all_data <- rbind(all_data, output)
#   }
#   
#   # Write processed data to CSV file
#   write.csv(all_data, file_name, row.names = FALSE)
# }


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
  filemap <- list('stl-target-rotation'='pointTargetSTL', 'stl-arc-rotation'='arcTargetSTL')
  identifiers <- getidentifiers()
  
  for (maxrot in c(45,60,90)) {
    
    for (condition in triallabels) {
      
      condition_data <- NA
      
      for (identifier in identifiers) {
        processed_data <- Preprocess_file(identifier, condition)
        processed_data$target <- list('stl-target-rotation'='point', 'stl-arc-rotation'='arc')[[condition]]
        if (max(processed_data$rotation) != maxrot) {
          next
        }
        if (is.data.frame(condition_data)) {
          condition_data <- rbind(condition_data, processed_data)
        } else {
          condition_data <- processed_data
        }
      }
      
      
      write.csv( x= condition_data,
                 file= sprintf('data/%s_%d.csv', filemap[[condition]], maxrot),
                 quote=FALSE,
                 row.names =FALSE) 
    }
  }

}

# # STL Prediction Function
# STLPredict <- function(par, rotations) {
#   
#   return((rotations * par['s']) * (dnorm(rotations, mean=0, sd= par['w']) / dnorm(0, mean=0, sd= par['w'])))
# }
# 
# # # Function to establish & return predicted adaptions for given fit data and rotations  
# # STLPredictions <- function(fits_data, rotations) {
# #   results <- data.frame(participant = character(), rotations = numeric(), predictions = numeric(), stringsAsFactors = FALSE)
# #   
# #   fits_data$slope <- as.numeric(fits_data$slope)
# #   fits_data$width <- as.numeric(fits_data$width)
# #   
# #   for (i in 1:nrow(fits_data)) {
# #     participant <- fits_data$participant[i]
# #     par <- c(s = fits_data$slope[i], w = fits_data$width[i])
# #     predictions <- STLPredict(par = par, rotations = rotations)
# #     
# #     participant_results <- data.frame(participant = rep(participant, length(rotations)),
# #                                       rotations = rotations,
# #                                       predictions = predictions)
# #     
# #     results <- rbind(results, participant_results)
# #   }
# #   
# #   return(results)
# # }    
# 
# 
# # STL Error Function
# STLErrors <- function(par, rotations, deviations) {
# 
#   # Generate model predictions using STLpredict function
#   model_predictions <- STLPredict(par, rotations)
# 
#   # Calculate squared errors between model predictions and actual deviations
#   squared_errors <- (deviations - model_predictions)^2
# 
#   # Compute the mean squared error (MSE)
#   MSE <- mean(squared_errors)
# 
#   return(MSE)
# }
# 
# # STL Grid Search Function
# STLGridsearch <- function(rotations, deviations) {
#   # Define parameter grids
#   s_values <- seq(0.1, 1.5, length.out = 15)  # Example range for parameter s
#   w_values <- seq(0, 60, length.out = 60)  # Example range for parameter w
# 
#   # Generate all combinations of parameters
#   parameter_combinations <- expand.grid(s = s_values, w = w_values)
# 
#   MSE <- apply(parameter_combinations, FUN=STLErrors, MARGIN=c(1), rotations = rotations, deviations = deviations)
# 
#   top_parameters <- parameter_combinations[order(MSE)[1:10], ]
# 
#   return(top_parameters)
# }
# 
# 
# STLOptimization <- function(top_parameters, rotations, deviations) { # run optimx on the best starting positions:
#   allfits <- do.call("rbind",
#                      apply( top_parameters,
#                             MARGIN=c(1),
#                             FUN=optimx::optimx,
#                             fn=STLErrors,
#                             method = "L-BFGS-B",
#                             lower = c(0.1, 1),
#                             upper = c(1.5, 60),
#                             rotations = rotations,
#                             deviations = deviations
#                      ))
# 
#   win <- allfits[order(allfits$value)[1], ]
#   return(unlist(win[1:2]))
# }

# # Returns individual STL Fit parameters for a list of participants
# STLIndividualFits <- function(identifiers_list, data_file = NULL) {
#   # Create an empty data frame to store the best fits for all participants
#   all_fits_df <- data.frame(participant = character(),
#                             slope = numeric(),
#                             width = numeric(),
#                             stringsAsFactors = FALSE)
#   
#   # Check if a data_file argument is provided
#   if (!is.null(data_file)) {
#     # Read the data file
#     data <- read.csv(data_file, stringsAsFactors = FALSE)
#     
#     # Loop through each participant
#     for (identifier in identifiers_list) {
#       # Extract data for the current participant from the data file
#       participant_data <- data[data$Participant_ID == identifier, ]
#       
#       # Extract rotations and deviations
#       rotations <- participant_data$Rotation
#       deviations <- participant_data$Response
#       
#       # Run grid search to find best parameters
#       top_parameters <- STLGridsearch(rotations, deviations)
#       
#       # Run optimization on the best starting positions
#       winning_par <- STLOptimization(top_parameters, rotations, deviations)
#       
#       # Add the participant id, slope, and width to the data frame
#       all_fits_df[nrow(all_fits_df) + 1, ] <- c(identifier, winning_par[1], winning_par[2])
#     }
#   } else {
#     # Loop through each participant
#     for (identifier in identifiers_list) {
#       # Preprocess data for the current participant
#       processed_data <- Preprocess_file(identifier)
#       
#       # Extract rotations and deviations
#       rotations <- processed_data$output$rotation
#       deviations <- processed_data$output$response
#       
#       # Run grid search to find best parameters
#       top_parameters <- STLGridsearch(rotations, deviations)
#       
#       # Run optimization on the best starting positions
#       winning_par <- STLOptimization(top_parameters, rotations, deviations)
#       
#       # Add the participant id, slope, and width to the data frame
#       all_fits_df[nrow(all_fits_df) + 1, ] <- c(identifier, winning_par[1], winning_par[2])
#     }
#   }
#   
#   # Set column names
#   colnames(all_fits_df) <- c("participant", "slope", "width")
#   
#   return(all_fits_df)
# }
# 
# 
# # Returns a single (or weighted) set of STL Fit parameters for a set of participant data
# STLFit <- function(participant, data_file) {
#   # Read raw data from CSV
#   data <- read.csv(data_file)
#   
#   # If participant is a single ID, convert it to a list
#   if (!is.vector(participant)) {
#     participant <- vector(participant)
#   }
#   
#   # Extract data for all participants
#   all_data <- NA
#   for (id in participant) {
#     if (is.data.frame(all_data)) {
#       all_data <- rbind(all_data, data[which(data$Participant_ID == id),])
#     } else {
#       all_data <- data[which(data$Participant_ID == id),]
#     }
#   }
#   
#   # Extract rotations and deviations from combined data
#   rotations <- all_data$Rotation
#   deviations <- all_data$Response
#   
#   # Run grid search to find best parameters
#   top_parameters <- STLGridsearch(rotations, deviations)
#   
#   # Run optimization on the best parameters
#   winning_par <- STLOptimization(top_parameters, rotations, deviations)
#   
#   # Create a data frame for the results
#   results_df <- data.frame(Participant_ID = "Combined",
#                            Slope = winning_par[1],
#                            Width = winning_par[2],
#                            stringsAsFactors = FALSE)
#   
#   return(results_df)
# }
# 
# # Baselines and derives Exponential Model Fits from participant Fixed Rotation data 
# STLExponentialFits <- function(ids) {
#   if (!is.vector(ids)) {
#     ids <- c(ids)  # Ensure ids is a vector
#   }
#   
#   results <- data.frame(participant = character(), lambda = numeric(), N0 = numeric(), stringsAsFactors = FALSE)
#   
#   for (id in ids) {
#     filename <- sprintf('data/%s_performance.csv', id)
#     df <- read.csv(filename, stringsAsFactors = FALSE)
#     
#     # Extract baseline data
#     bias <- median(df$reachdeviation_deg[which(df$label == 'fixed-rotation-baseline')[c(11:60)]])
#     learning <- df$reachdeviation_deg[which(df$label == 'fixed-rotation')] - bias
#     
#     # Derive Exponential Fits
#     fit <- Reach::exponentialFit(signal = learning)
#     
#     # Initialize a data frame to store Exponential Fit data
#     participant_result <- data.frame(participant = id, lambda = fit['lambda'], N0 = fit['N0'])
#     
#     # Store Exponential Fit Data
#     results <- rbind(results, participant_result)
#   }
#   
#   return(results)
# }
#     
# 
# # Obtains graph-able Exponential Model values for participants
# STLExponentialModel <- function(participants_results, mode = 'learning', setN0 = NULL, points) {
#     
#     participants_results$lambda <- as.numeric(participants_results$lambda)
#     participants_results$N0 <- as.numeric(participants_results$N0)
#     
#     # Initialize an empty data frame to hold all the plot data
#     all_plot_data <- data.frame()
#     
#     # Iterate over each participant's results
#     for (i in 1:nrow(participants_results)) {
#       participant_result <- participants_results[i, ]
#       lambda <- participant_result$lambda
#       N0 <- participant_result$N0
#       participant <- participant_result$participant
#       
#       # Generate the time points (for example, 0 to 100)
#       time_points = seq(0, points, by = 1)
#       
#       # Generate the model values using the exponentialModel function
#       model_output <- Reach::exponentialModel(par = c('lambda' = lambda, 'N0' = N0), timepoints = time_points, mode = mode, setN0 = setN0)
#       
#       # Create a data frame for the current participant's plot data
#       plot_data <- data.frame(time = model_output$trial, value = model_output$output, participant = participant)
#       
#       # Combine with the overall plot data
#       all_plot_data <- rbind(all_plot_data, plot_data)
#     }
#     return(all_plot_data)
# }
#     
# # Function to plot a participant's observed STL reach data along with average & predicted deviation lines
# PlotReachdata <- function(identifiers_list, filename, predictions_data) {
#   
#   dir_path <- "~/Desktop/Masters/Thesis/PredictiveSTL/Plots/Raw Data Plots/"
#   
#   pdf_name <- paste0(dir_path, filename, "_Reach_Data.pdf")
#   
#   if (!dir.exists(dir_path)) {
#     dir.create(dir_path, recursive = TRUE)
#   }
#   
#   pdf(pdf_name, width = 8.5, height = 11)  # Open PDF device
#   
#   layout(mat=matrix(c(1:6),byrow=TRUE,nrow=3))
#   
#   data_file <- read.csv("~/Desktop/Masters/Thesis/PredictiveSTL/Data Files/raw_data.csv", stringsAsFactors = FALSE)
#   
#   # Loop through each participant
#   for (identifier in identifiers_list) {
#     # Extract data for the current participant from the data file
#     participant_data <- data_file[data_file$Participant_ID == identifier, ]
#     
#     # Extract rotations and deviations
#     rotations <- participant_data$Rotation
#     deviations <- participant_data$Response
#     
#     # Calculate Mean Deviation per rotation value
#     mean_deviation <- participant_data %>%
#       group_by(Rotation) %>%
#       summarize(mean_deviation = mean(Response))
#     
#     # Extract predicted data for the current participant
#     participant_predictions <- predictions_data %>%
#       filter(participant == identifier)
#     
#     p <- ggplot(participant_data, aes(x = Rotation, y = Response)) +
#       geom_point(color = "black") +
#       
#       # Plot mean data points with red line
#       geom_line(data = mean_deviation, aes(x = Rotation, y = mean_deviation), color = "red", alpha = 0.6) +
#       
#       # Plot predicted data points with blue line
#       geom_line(data = participant_predictions, aes(x = rotations, y = predictions), color = "blue", alpha = 0.6) +
#       
#       # Add a dotted grey line at y = 0
#       geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +    
#       
#       labs(title = identifier,
#            x = "Rotation in Degrees(°)", 
#            y = "Deviation in Degrees(°)") +
#       theme_minimal() +
#       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#     
#     # Print the plot to the PDF device
#     print(p)
#    
#   }
#   
#   dev.off()  # Close PDF device
# 
# }
# 
# # Function to plot exponential model values 
# STLExponentialPlot <- function(data, rot, vert_line = NULL, trials = NULL){
#   
#   main_title <- paste(paste0(rot,"°"), "Adaptation Over", trials, "Trials")
#   
#   # Allows to plot for a certain number of trials
#   if (is.null(trials)){
#   
#   ggplot(data, aes(x = time, y = value, color = participant)) +
#     geom_line(alpha = 0.7) +
#       labs(title = main_title,
#            x = "Trials",
#            y = "Model Value") +
#       theme_minimal() +
#       guides(color = "none") +
#       geom_vline(xintercept = vert_line, linetype = "dotted", color = "black")
#     
#     
#   } else{
#     # Add the trials column, resetting at 0 for each new participant
#       data <- data %>%
#         group_by(participant) %>%
#         mutate(trials = row_number() - 1) %>%
#         ungroup()
#       
#       # Filter data to include only up to the specified number of trials
#       data <- data %>%
#         filter(trials <= trials)
#       
#       ggplot(data, aes(x = trials, y = value, color = participant)) +
#         geom_line(alpha = 0.7) +
#         labs(title = main_title,
#              x = "Trials",
#              y = "Model Value") +
#         theme_minimal() +
#         guides(color = "none") +
#         geom_vline(xintercept = vert_line, linetype = "dotted", color = "black") +
#         scale_x_continuous(limits = c(0, trials))
#       
#   }
# }
# 
# # Function to plot predicted STL adaptation values
# STLPredictPlot <- function(data, vert_line = NULL){
#   
# main <- "Predicted STL Adaption"
#   
# ggplot(data, aes(x = rotations, y = predictions, group = participant)) +
#   geom_line(color = "red", alpha = 0.6) +  # Lower alpha for more transparency
#   labs(title = main,
#        x = "Rotations in Degrees(°)",
#        y = "Predicted Adaptation in Degrees(°)") +
#   theme_minimal() +
#   geom_vline(xintercept = vert_line, linetype = "dotted", color = "black") 
# }   

# load processed data ----

loadSTLdata <- function(targets=c('arc','point'),
                        average=NULL,
                        maxrots=c(45,60),
                        superblocks=NULL) {
  
  if (is.null(superblocks)) {
    careaboutsuperblocks <- FALSE
    superblocks <- c(1:5)
  } else {
    careaboutsuperblocks <- TRUE
  }
  df <- NA
  for (maxrot in maxrots) {
    for (target in targets) {
      tdf <- read.csv(sprintf('data/%sTargetSTL_%d.csv', target, maxrot), stringsAsFactors = FALSE)
      tdf <- tdf[which(tdf$superblock %in% superblocks),]
      tdf$maxrotation <- maxrot
      if (is.data.frame(df)) {
        df <- rbind(df,tdf)
      } else {
        df <- tdf
      }
    }
  }

  if (!(is.null(average))) {
    if (careaboutsuperblocks) {
      df <- aggregate(response ~ target + participant + rotation + maxrotation + superblock, data=df, FUN=average)
    } else {
      df <- aggregate(response ~ target + participant + rotation + maxrotation, data=df, FUN=average)
    }
  }
  
  return(df)

}

getMaxRot <- function(df) {
  
  idx <- which(df$label == 'stl-target-rotation')
  
  return(max(df$rotation[idx]))
  
}

extractLongExpData <- function(df) {
  
  baseline_idx <- which(df$label == 'fixed-rotation-baseline')
  exposure_idx <- which(df$label == 'fixed-rotation')
  aftereffects_idx <- which(df$label == 'fixed-rotation-aftereffects')
  
  # print(c(length(baseline_idx),length(exposure_idx),length(aftereffects_idx)))
  
  baseline <- df$reachdeviation_deg[baseline_idx[31:60]]
  baseline <- baseline[which(abs(baseline) < 60)]
  baseline <- mean(baseline)
  # print(baseline)
  
  exposure <- df$reachdeviation_deg[exposure_idx] - baseline
  exposure[which(abs(exposure) > 80)] <- NA
  
  aftereffects <- df$reachdeviation_deg[aftereffects_idx] - baseline
  aftereffects[which(abs(aftereffects) > 80)] <- NA
  
  return( list( 'baseline'      = baseline,
                'exposure'      = exposure,
                'aftereffects'  = aftereffects) )
  
}

getLongTimeCourses <- function(maxrots = c(45,60,90)) {
  
  IDs <- getidentifiers()
  
  dfs <- list('45' = NA,
              '60' = NA,
              '90' = NA)
  
  for (IDno in c(1:length(IDs))) {
    
    ID <- IDs[IDno]
    
    pdf  <- read.csv(sprintf('data/%s_performance.csv',ID), stringsAsFactors = FALSE)
    mrot <- getMaxRot(pdf)
    ldf  <- extractLongExpData(pdf)
    
    pdf <- data.frame( participant = rep(ID,160),
                       trial       = c(1:160),
                       reachdeviation_deg = ldf[['exposure']],
                       maxrot      = rep(mrot,160) )
    
    idx <- sprintf('%d', mrot)
    
    if (is.data.frame(dfs[[idx]])) {
      dfs[[idx]] <- rbind(dfs[[idx]], pdf)
    } else {
      dfs[[idx]] <- pdf
    }
    
  }
  
  return(dfs)
  
}

getFirstTrialData <- function() {
  
  # we need the exponential fit data:
  exponentials <- read.csv('data/exponential_fits.csv', stringsAsFactors = FALSE)
  exponentials <- exponentials[which(exponentials$phase == 'learning'),]
  exponentials$prediction <- exponentials$lambda * exponentials$N0
  
  # and the long time course data itself:
  exposure <- getLongTimeCourses()
  exposure <- rbind(exposure[['45']], exposure[['60']], exposure[['90']])
  # exposure <- exposure[,-which(names(exposure) %in% c('maxrot'))]
  
  exposure <- exposure[which(exposure$trial == 2),]
  
  firstTrials <- merge(exponentials, exposure, by=c('participant', 'maxrot'))
  
  return(firstTrials)
  
}