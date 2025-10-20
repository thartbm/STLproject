
cappedModel <- function(par, rotations) {
  
  slope <- par['r']
  cap <- par['c']
  
  return(pmin(cap, slope * rotations))
  
}

cappedMSE <- function(par, rotations, deviations) {
  
  predictions <- cappedModel(par, rotations)
  
  squared_errors <- (predictions - deviations)^2
  
  MSE <- mean(squared_errors)
  
  return(MSE)
  
}

cappedGridSearch <- function(rotations, deviations) {

  # r_values <- seq(0.1, 4.9, length.out = 25)  # Example range for parameter r
  
  r_values <- seq(0.05, .95, length.out = 25)  # Example range for parameter r
  c_values <- seq(1, 25, length.out = 25)  # Example range for parameter w
    
  # Generate all combinations of parameters
  parameter_combinations <- expand.grid(r = r_values, c = c_values)
    
  MSE <- apply(parameter_combinations, FUN=cappedMSE, MARGIN=c(1), rotations = rotations, deviations = deviations)
    
  top_parameters <- parameter_combinations[order(MSE)[1:10], ]
    
  return(top_parameters)
  
}

cappedOptimization <- function(top_parameters, rotations, deviations) { # run optimx on the best starting positions:
  allfits <- do.call("rbind",
                     apply( top_parameters,
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=cappedMSE,
                            method = "L-BFGS-B", 
                            lower = c(0.001, 0), 
                            # upper = c(5.000, 25),
                            upper = c(1.000, 25),
                            rotations = rotations,
                            deviations = deviations
                     ))
  
  win <- allfits[order(allfits$value)[1], ]
  return(unlist(win[1:2]))
}

cappedFit <- function(participants, data) {

  if (!is.vector(participants)) {
    participants <- vector(participants)
  }
  
  all_data <- NA
  for (id in participants) {
    if (is.data.frame(all_data)) {
      all_data <- rbind(all_data, data[which(data$participant == id),])
    } else {
      all_data <- data[which(data$participant == id),]
    }
  }
  
  rotations <- all_data$rotation
  deviations <- all_data$response
  
  top_parameters <- cappedGridSearch(rotations, deviations)
  
  winning_par <- cappedOptimization(top_parameters, rotations, deviations)
  
  return(winning_par)
}