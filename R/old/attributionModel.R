library(Reach)
library(optimx)

# STL Prediction Function
STLPredict <- function(par, rotations) {

    return((rotations * par['s']) * (dnorm(rotations, mean=0, sd= par['w']) / dnorm(0, mean=0, sd= par['w'])))
  
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
  s_values <- seq(0.0625, 1.0, length.out = 16)  # Example range for parameter s
  w_values <- seq(1, 40, length.out = 18)  # Example range for parameter w # lower value CAN NOT BE 0
  
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
                            lower = c(0.001, 0.1),
                            upper = c(1.000, 40),
                            rotations = rotations,
                            deviations = deviations
                     ))
  
  win <- allfits[order(allfits$value)[1], ]
  return(unlist(win[1:2]))
}

# Returns a single (or weighted) set of STL Fit parameters for a set of participant data
STLFit <- function(participants, data) {

  # If participant is a single ID, convert it to a list
  if (!is.vector(participants)) {
    participants <- vector(participants)
  }

  # Extract data for all participants
  all_data <- NA
  for (id in participants) {
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
  return( STLOptimization(top_parameters, rotations, deviations) )

}