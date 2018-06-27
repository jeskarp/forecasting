# R script for forecasting an epidemic's trajectory

#########################################
## Functions for input into prediction ##
#########################################

# A discretised serial interval
si_disc <- function(mu, cv) {
  library(epitrix)
  sim_si <- gamma_mucv2shapescale(mu, cv)
  si <- distcrete("gamma", shape = sim_si$shape, scale = sim_si$shape, w = 0, interval = 1)
  return(si)
}

# A simulated outbreak
sim_outbreak <- function(si, obs_time, R0) {
  # Simulate outbreak
  # set.seed(1)
  sim_test <- simOutbreak(R0 = R0, infec.curve = si, n.hosts = 1000000, duration = obs_time, seq.length = 10,
                          stop.once.cleared = FALSE)
  # Put the output that I care about into a data.frame
  sim_outbreak <- data.frame(id = sim_test$id,
                             inf_id = sim_test$ances,
                             onset = sim_test$onset)
  return(sim_outbreak)
}

# Aggregated data required for running projections
full_data <- function() {
  # Serial interval
  mu <- 19 # SI mean
  # cv <- 0.75 # SI coefficient of variance
  sigma <- 11.0 # SI standard deviation 
  
  cv <- sigma / mu # SI coefficient of variance
  # sigma <- cv * mu # SI standard deviation
  
  serial_interval <- si_disc(mu, cv)
  
  # Number of days you'd like to forecast for
  # make calibration + projection 8 delta - where delta is SI mu/2
  delta <- round(mu/2) # smallest unit of time over which I calibrate/project
  no_chunks <- 8 # number of data chunks that will be split by delta
  
  # Simulated outbreak
  obs_time <- delta * no_chunks # how long simulation is run for
  R0 <- 1.71 # R0 of the outbreak
  outbreak_data <- sim_outbreak(serial_interval$d(0:30), obs_time, R0)
  
  # Number of trajectories for the projection
  n_traj <- 10000

  full_data <- list("outbreak_data" = outbreak_data, "serial_interval" = serial_interval, 
                    "mu" = mu, "sigma" = sigma, "delta" = delta, "no_chunks" = no_chunks, 
                    "n_traj" = n_traj)
  
  return(full_data)
}

# Listing the calibration and projection combinations
split_data <- function(outbreak_data) {
  all_combos <- expand.grid(1:(outbreak_data$no_chunks - 1), 8)
  all_combos <- setNames(all_combos, c("calibration","projection"))
  return(all_combos)
}

##############################
## Functions for projection ##
##############################

# Make daily incidence of the outbreak data
observed_incidence <- function(linelist) {
  library(incidence)
  observed_incidence <- incidence(linelist$onset, interval = 1)
  return(observed_incidence)
}

R_estimate <- function(observed_incidence, serial_interval) {
  # Estimate R
  library(earlyR)
  # The R estimation
  R <- get_R(observed_incidence, si = serial_interval, 
             max_R = 10)
  return(R)
}

projection <- function(full_data, obs_incidence, cutoff_time, proj_start, proj_end) {
  # Forecast incidence
  library(projections)

  # number of days that I hid data for i.e. number of days I'm projecting for
  hidden_days <- proj_end - cutoff_time
  
  # Call the R_estimate function
  R <- R_estimate(obs_incidence[1:cutoff_time, ], full_data$serial_interval)
  
  # Project
  proj <- project(obs_incidence[1:cutoff_time, ], R = sample_R(R, 1000), si = full_data$serial_interval, 
                  n_sim = full_data$n_traj, n_days = hidden_days, R_fix_within = TRUE)
  return(proj)
}

##########################
## Multiple simulations ##
##########################

# Make multiple simulations and project for them
# multi_simulation <- function(n_sim) {
#   for (i in 1:n_sim) {
      # Will need to specify where in the dataframe this particular output will be saved
#     simulation <- output()
#   }
#   return(simulation)
# }

########################
## Prediction metrics ##
########################

## Reliability

# need the true datapoints for a given time and the predictions for the same time

# The cumulative probability distribution for time t: F_t
cumulative_poisson <- function(data, pred){
  ppois(data, mean(pred))
}

reliability <- function(data, pred) {
  library(incidence)
    reliability <- array(NA, dim = c(nrow(data))) # nrow(data) is how many days we have hidden data for
    for (i in 1:nrow(data)){
      reliability[i] <- cumulative_poisson(data[i], pred[i, ])
    }
    return(reliability)
  # anderson-darling <- function(data, pred){
  #   anderson-darling <- array(NA, dim = c(hidden_days))
  #   for (i in 1:hidden_days){
  #   ad.test()
  #   }
  #   return(anderson-darling)
  # }
}
# model_calibration_sim_1_all <- reliability(sim_hidden_incidence_1_all$counts, proj_1[1:(serial_int * 3), ])



## Sharpness

# Function for calculating forecasted incidence median
forecast_median <- function(x){
  return(median(x))
}

# Function for calculating the standard deviations
# Want it to return MADM(y)
MADM <- function(x){
  traj_diff <- array(NA, dim = c(length(x))) # array for storing SDs
  # For each trajectory for this day, calculate the SD
  for (i in 1:length(x)){
    # save the absolute difference for each trajectory
    traj_diff[i] <-  abs(x[i] - forecast_median(x))
  }
  MADM <- median(as.vector(traj_diff))
  return(MADM)
}

# Function for calculating the normalised median absolute deviation
# Want it to return St(Ft)
sharpness <- function(x){
  # array for storing daily sharpness
  sharpness <- array(NA, dim = c(nrow(x)))
  for (i in 1:nrow(x)){
    sharpness[i] <- 1 - (MADM(as.vector(x[i, ])) / forecast_median(as.vector(x[i, ]))) 
  }
  return(sharpness)
}

# model_sharpness_sim_1_all <- sharpness(proj_1[1:(serial_int * 3), ])



## Bias

# Calculate the difference between datapoints
# Feed in the true datapoint and predictions for a given time t
difference <- function(data, pred){
  # Calculate difference between data point and sample
  difference <- array(NA, dim = c(length(pred)))
  for (i in 1:length(pred)){
    difference[i] <- pred[i] - data
  }
  return(difference)
}

# Do the Heaviside function with half-maximum convention
# Feed in the differences between the predictions and true values for a given time t
heaviside <- function(difference) {
  heaviside <- array(NA, dim = c(length(difference)))
  for (i in 1:length(difference)) {
    if (is.na(difference[i]) == TRUE) {
      heaviside[i] <- NA
    } else if (difference[i] > 0) {
      heaviside[i] <- 1
    } else if (difference[i] < 0) {
      heaviside[i] <- 0
    } else {
      heaviside[i] <- 0.5
    }
  }
  return(heaviside)
}

# Calculate the expected prediction value for a given time t
expected <- function(h_values, pred){
  mean_heaviside <- mean(h_values, na.rm = TRUE)
  return(mean_heaviside)
}

# Calculate bias
bias <- function(data, pred){
  bias <- array(NA, dim = c(nrow(data)))
  for (i in 1:nrow(data)){
    bias[i] <- 2 * (expected(heaviside(difference(data[i], as.vector(pred[i, ]))), as.vector(pred)) - 0.5)
  }
  return(bias)
}

# Get single bias score
bias_score <- function(data, pred){
  daily_bias <- bias(data, pred)
  bias_score <- sum(daily_bias) / nrow(data)
  return(bias_score)
}
# model_bias_sim_1_all <- bias(sim_hidden_incidence_1_all$counts, proj_1[1:(serial_int * 3), ])
# bias_score(sim_hidden_incidence_1_all$counts, proj_1[1:(serial_int * 3), ])



## Root-mean-square error

# Calculate squared difference between true data and prediction for a given timepoint
squared_diff <- function(data, pred){
  squared_diff <- array(NA, dim = c(length(pred)))
  for (i in 1:length(pred)){
    squared_diff[i] <- (data - pred[i])^2
  }
  return(squared_diff)
}

# Take the mean of the squared differences
mean_squared_diff <- function(data, pred){
  mean_squared_diff <- sum(squared_diff(data, pred)) / length(pred)
  return(mean_squared_diff)
}

# Calculates RMSE for each forecasted timepoint
rmse <- function(data, pred){
  rmse <- array(NA, dim = c(length(data)))
  for (i in 1:length(data)){
    rmse[i] <- sqrt(mean_squared_diff(data[i], as.vector(pred[i, ])))
  }
  return(rmse)
}
# model_rmse_sim_1_all <- rmse(sim_hidden_incidence_1_all$counts, proj_1[1:(serial_int * 3), ])
# sum_rmse_1 <- sum(model_rmse_sim_1_all)

################################
## Functions for data storage ##
################################

hash_function <- function(projected_incidence) {
  library("digest")
  hash <- digest(projected_incidence)
  return(hash)
}

chunk_projection <- function(full_data, combo_list) {
  
  for (i in 1:nrow(combo_list)) {
    cutoff_time <- full_data$delta * combo_list[i, 1] # the time at which observed data stops
    proj_end <- full_data$delta * combo_list[i, 2] # the time at which projection ends
    proj_start <- cutoff_time + 1 # the time at which projection starts
    proj_window <- projection(full_data, obs_incidence, cutoff_time, proj_start, proj_end)
  }
}

#################################################
## The projection and prediction metric output ##
#################################################

output <- function(n_sim) {
  library(dplyr)
  
  # Set seed for reproducibility
  # set.seed(1)
  sim_gen_data <- full_data() # has all the original data I need
  
  # Create a directory into which I store this simulation
  sim_hash <- digest(sim_gen_data) # unique identifier for simulation
  dir.create(paste("/home/evelina/Development/forecasting/simulations/", sim_hash, sep = ""))
  
  # Set this directory as the working directory so files end up in the right place
  setwd(paste("/home/evelina/Development/forecasting/simulations/", sim_hash, sep = ""))
  
  # Projections and prediction metrics
  obs_incidence <- observed_incidence(sim_gen_data$outbreak_data) # has the incidence object for the outbreak
  combo_list <- split_data(sim_gen_data) # all the combinations of deltas that I want to project for
  
  for (i in 1:nrow(combo_list)) {
    cutoff_time <- sim_gen_data$delta * combo_list[i, 1] # the time at which observed data stops
    proj_end <- sim_gen_data$delta * combo_list[i, 2] # the time at which projection ends
    proj_start <- cutoff_time + 1 # the time at which projection starts
    proj_window <- projection(sim_gen_data, obs_incidence, cutoff_time, proj_start, proj_end) # projection for the time window
    
    # plot(obs_incidence) %>% add_projections(proj_window[1:10])
    
    # Calculate prediction metrics
    proj_rel <- reliability(obs_incidence[proj_start:proj_end, ]$counts, proj_window)
    proj_sharp <- sharpness(proj_window)
    proj_bias <- bias(obs_incidence[proj_start:proj_end, ]$counts, proj_window)
    proj_bias_score <- bias_score(obs_incidence[proj_start:proj_end, ]$counts, proj_window)
    proj_rmse <- rmse(obs_incidence[proj_start:proj_end, ]$counts, proj_window)
    proj_rmse_sum <- sum(proj_rmse, na.rm = TRUE)
    
    # make an array for storing the prediction metrics for a given projection
    proj_metrics <- data.frame(dataset = sim_hash,
                               cali_window_size = cutoff_time,
                               proj_window_size = sim_gen_data$delta,
                               disease = "ebola",
                               days_since_data = c(1:(proj_end - cutoff_time)),
                               reliability = proj_rel,
                               sharpness = proj_sharp,
                               bias = proj_bias,
                               bias_score = proj_bias_score,
                               rmse = proj_rmse,
                               sum_rmse = proj_rmse_sum)
    
    # save this projection for this window
    save(proj_window, file = paste("proj_window_", i, ".RData", sep = ""))
    
    # rbind proj_metrics to a bigger data.frame that will hold all the projections' relevant outputs for a simulation
    if (i == 1) {
      full_proj_metrics <- proj_metrics 
    } else {
      full_proj_metrics <- dplyr::bind_rows(full_proj_metrics, proj_metrics)
    }
  }
  
  # Save important information
  write.csv(full_proj_metrics, file = "full_proj_metrics.csv")
  save(sim_gen_data, file = "sim_generating_data.RData")
  
  return(paste("Finished projection", n_sim, sep = " "))
}

multi_output <- function(n_simulations) {
  for (i in 1:n_simulations){
    output(i)
  }
  return("Finished all projections")
}
