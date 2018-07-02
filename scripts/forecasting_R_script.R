# R script for forecasting an epidemic's trajectory

#########################################
## Functions for input into prediction ##
#########################################

# A discretised serial interval
si_disc <- function(mu, cv) {
  library(epitrix)
  library(distcrete)
  sim_si <- gamma_mucv2shapescale(mu, cv)
  si <- distcrete("gamma", shape = sim_si$shape, scale = sim_si$shape, w = 0, interval = 1)
  return(si)
}

# A simulated outbreak
sim_outbreak <- function(si, obs_time, R0) {
  library(outbreaker)
  # Simulate outbreak
  set.seed(3)
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
  n_traj <- 10 # 10000

  full_data <- list("outbreak_data" = outbreak_data, "serial_interval" = serial_interval, 
                    "mu" = mu, "sigma" = sigma, "delta" = delta, "no_chunks" = no_chunks, 
                    "n_traj" = n_traj)
  
  return(full_data)
}

# Listing the calibration and projection combinations
split_data <- function(outbreak_data) {
  # project for maximum 4 deltas
  # latest projection can be last chunk of no_chunks
  all_combos <- expand.grid(1:(outbreak_data$no_chunks - 1), outbreak_data$no_chunks)
  all_combos <- setNames(all_combos, c("calibration","projection"))
  return(all_combos)
}

##############################
## Functions for projection ##
##############################

# Make daily incidence of the outbreak data
observed_incidence <- function(linelist, delta, no_chunks) {
  library(incidence)
  observed_incidence <- incidence(linelist$onset, interval = 1, last_date = (delta * no_chunks))
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

######################################
## Functions for prediction metrics ##
######################################

# x_t - true data for a given outbreak

## Reliability

# need the true datapoints for a given time and the predictions for the same time

# The cumulative probability distribution for time t: F_t
cumulative_poisson <- function(x_t, pred){
  ppois(x_t, mean(pred))
}

reliability <- function(x_t, pred) {
  library(incidence)
    reliability <- array(NA, dim = c(nrow(x_t))) # nrow(data) is how many days we have hidden data for
    for (i in 1:nrow(x_t)){
      reliability[i] <- cumulative_poisson(x_t[i], pred[i, ])
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



## Sharpness

# Function for calculating forecasted incidence median
# forecast_median <- function(x_t){
#   return(median(x_t))
# }

# Function for calculating the standard deviations
# Want it to return MADM(y)
MADM <- function(pred){
  traj_diff <- array(NA, dim = c(length(pred))) # array for storing SDs
  # For each trajectory for this day, calculate the SD
  for (i in 1:length(pred)){
    # save the absolute difference for each trajectory
    traj_diff[i] <-  abs(pred[i] - median(pred)) # median is the forecast median
  }
  MADM <- median(as.vector(traj_diff))
  return(MADM)
}

# Function for calculating the normalised median absolute deviation
# Want it to return St(Ft)
sharpness <- function(pred){
  # array for storing daily sharpness
  sharpness <- array(NA, dim = c(nrow(pred)))
  for (i in 1:nrow(pred)){
    sharpness[i] <- 1 - (MADM(as.vector(pred[i, ])) / median(as.vector(pred[i, ]))) # median is the forecast median
  }
  return(sharpness)
}



## Bias

# Calculate the difference between datapoints
# Feed in the true datapoint and predictions for a given time t
difference <- function(x_t, pred){
  # Calculate difference between data point and sample
  difference <- array(NA, dim = c(length(pred)))
  for (i in 1:length(pred)){
    difference[i] <- pred[i] - x_t
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
bias <- function(x_t, pred){
  bias <- array(NA, dim = c(nrow(x_t)))
  for (i in 1:nrow(x_t)){
    bias[i] <- 2 * (expected(heaviside(difference(x_t[i], as.vector(pred[i, ]))), as.vector(pred)) - 0.5)
  }
  return(bias)
}

# Get single bias score
bias_score <- function(x_t, pred){
  daily_bias <- bias(x_t, pred)
  bias_score <- sum(daily_bias) / nrow(x_t)
  return(bias_score)
}



## Root-mean-square error

# Calculate squared difference between true data and prediction for a given timepoint
squared_diff <- function(x_t, pred){
  squared_diff <- array(NA, dim = c(length(pred)))
  for (i in 1:length(pred)){
    squared_diff[i] <- (x_t - pred[i])^2
  }
  return(squared_diff)
}

# Take the mean of the squared differences
mean_squared_diff <- function(x_t, pred){
  mean_squared_diff <- sum(squared_diff(x_t, pred)) / length(pred)
  return(mean_squared_diff)
}

# Calculates RMSE for each forecasted timepoint
rmse <- function(x_t, pred){
  rmse <- array(NA, dim = c(length(x_t)))
  for (i in 1:length(x_t)){
    rmse[i] <- sqrt(mean_squared_diff(x_t[i], as.vector(pred[i, ])))
  }
  return(rmse)
}

#################################################
## The projection and prediction metric output ##
#################################################

output <- function(n_sim) {
  library(dplyr)
  library(incidence)
  library(digest)
  
  # Set seed for reproducibility
  # set.seed(1)
  sim_gen_data <- full_data() # has all the original data I need
  
  # Create a directory into which I store this simulation
  sim_hash <- digest(sim_gen_data) # unique identifier for simulation
  dir.create(paste("/home/evelina/Development/forecasting/simulations/", sim_hash, sep = ""))
  
  # Set this directory as the working directory so files end up in the right place
  setwd(paste("/home/evelina/Development/forecasting/simulations/", sim_hash, sep = ""))
  
  # Projections and prediction metrics
  obs_incidence <- observed_incidence(sim_gen_data$outbreak_data, sim_gen_data$delta, sim_gen_data$no_chunks) # has the incidence object for the outbreak
  combo_list <- split_data(sim_gen_data) # all the combinations of deltas that I want to project for
  
  # Save plot of incidence curve
  pdf("incidence.pdf", width = 7, height = 5)
    print(plot(obs_incidence))
  dev.off()    
  
  # A for loop for doing all the projections and calculating prediction metrics
  for (i in 1:nrow(combo_list)) {
    cutoff_time <- sim_gen_data$delta * combo_list[i, 1] # the time at which observed data stops
    proj_end <- sim_gen_data$delta * combo_list[i, 2] # the time at which projection ends
    proj_start <- cutoff_time + 1 # the time at which projection starts
    
    # if (i == 7) {
    #   print("We should have NAs here")
    # }
    
    proj_window <- projection(sim_gen_data, obs_incidence, cutoff_time, proj_start, proj_end) # projection for the time window
    
    # print(plot(obs_incidence) %>% add_projections(proj_window))
    
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
}

multi_output <- function(n_simulations) {
  for (i in 1:n_simulations){
    output(i)
  }
  return("Finished all projections")
}

output(1)