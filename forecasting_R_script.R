# R script for forecasting an epidemic's trajectory

####################
## Required input ##
####################

# Aggregated daily case data
data <- function() {
  daily_incidence <- c(1, 2, 3, 2, 3, 4, 5, 4, 7) # Actual data would go here
  length_obs <- 3 # How many days' worth of data you're observing for, the rest is hidden

  # Number of days you'd like to forecast for
  # Serial interval
  serial_interval <- c(1.2, 1.3, 2, 2.3, 2, 0, 0, 0, 0)
  mu <- 2 # This is where serial interval mean goes
  sigma <- 1 # this is where the serial interval SD goes 

  # The maximum value R can take for estimates
  highest_R <- 10

  # Number of trajectories for the projection
  n_traj <- 10 # 10000

  # Number of projections
  n_proj <- 10 # 1000
  
  data <- list("daily_incidence" = daily_incidence, "length_obs" = length_obs, "serial_interval" = serial_interval, "mu" = mu, "sigma" = sigma, "highest_R" = highest_R, "n_traj" = n_traj, "n_proj" = n_proj)
  
  return(data)
}
###################
## The functions ##
###################

observed_incidence <- function(daily_incidence, length_obs) {
  # Make incidence into an incidence object
  library(incidence)
  observed_incidence <- as.incidence(c(daily_incidence[1:length_obs])) # The data that you pretend you actually have
  return(observed_incidence)
}

R_estimate <- function(observed_incidence, mu, sigma, highest_R) {
  # Estimate R
  library(earlyR)
  # The R estimation
  R <- get_R(observed_incidence, si_mean = mu, si_sd = sigma, 
             max_R = highest_R)
  return(R)
}

projection <- function(daily_incidence, length_obs, serial_interval, mu, sigma, highest_R, n_traj) {
  # Forecast incidence
  library(projections)

  # number of days that I hid data for i.e. number of days I'm projecting for
  end <- length(daily_incidence) # when does the data end
  hidden_days <- end - length_obs
  
  # Call the observed_incidence and R_estimate functions
  observed_incidence <- observed_incidence(daily_incidence, length_obs)
  R <- R_estimate(observed_incidence, mu, sigma, highest_R)
  
  # Serial interval distribution already provided with the EpiEstim data
  proj <- project(observed_incidence, R = sample_R(R, 1000), si = serial_interval, 
                  n_sim = n_traj, n_days = hidden_days, R_fix_within = TRUE)
  return(proj)
}

# Do many projections
multi_projection <- function(data) {
  end <- length(data$daily_incidence) # when does the data end
  hidden_days <- end - data$length_obs
  
  multi_projection <- array(NA, c(hidden_days, data$n_traj, data$n_proj)) # create array for storing each projection
  
  for (i in 1:data$n_proj){
    single_projection <- projection(data$daily_incidence, data$length_obs, data$serial_interval, data$mu, data$sigma, data$highest_R, data$n_traj)
    multi_projection[, , i] <- single_projection   
  }
  return(multi_projection)
}

hash_function <- function(projected_incidence) {
  library("digest")
  hash <- digest(projected_incidence)
  return(hash)
}

reformat_dataset <- function(multi_projection, data) {
  dataset <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  for (i in 1:data$n_proj){
    for (j in 1:data$n_traj){
      hash <- array(hash_function(multi_projection[ , j, i]), c(length(multi_projection[ , j, i]), 1))
      time_seq <- seq(1, length(multi_projection[ , j, i]), 1)
      timepoint <- array(time_seq, c(length(multi_projection[ , j, i]), 1))
      hash_time <- cbind(hash, timepoint)
      trajectory_data <- cbind(hash_time, multi_projection[ , j, i])
      dataset <- rbind(dataset, trajectory_data)
    }
  }
  colnames(dataset) <- c("hash", "days_since_data", "num_cases")
  return(dataset)
}

################
## The output ##
################
output <- function() {
  # Set seed for reproducibility
  set.seed(1)
  data <- data()
  test <- multi_projection(data)
  reformat_test <- reformat_dataset(test, data)
  # write.csv(test, file = "test.csv")
  return(reformat_test)
}
# Daily case data for each simulation
# Days since data
# Hash
