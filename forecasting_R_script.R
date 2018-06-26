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
  set.seed(1)
  sim_test <- simOutbreak(R0 = R0, infec.curve = si, n.hosts = 1000000, duration = obs_time, seq.length = 10, stop.once.cleared = FALSE)
  # Put the output that I care about into a data.frame
  sim_outbreak <- data.frame(id = sim_test$id,
                             inf_id = sim_test$ances,
                             onset = sim_test$onset)
  return(sim_outbreak)
}

# Aggregated data required for running projections
full_data <- function() {
  # daily_incidence <- c(1, 2, 3, 2, 3, 4, 5, 4, 7) # Actual data would go here
  length_obs <- 3 # How many days' worth of data you're observing for, the rest is hidden

  # Serial interval
  mu <- 2 # SI mean
  cv <- 0.75 # SI coefficient of variance
  sigma <- cv * mu # SI standard deviation
  serial_interval <- si_disc(mu, cv)
  
  # Simulated outbreak
  obs_time <- 61 # how long simulation is run for
  R0 <- 1.5 # R0 of the outbreak
  outbreak_data <- sim_outbreak(serial_interval$d(0:30), obs_time, R0)
  
  # Number of days you'd like to forecast for
  # make calibration + projection 8 delta - where delta is SI mu/2
  delta <- round(mu/2) # smallest unit of time over which I calibrate/project
  no_chunks <- 8 # number of data chunks that will be split by delta
  
  # Number of trajectories for the projection
  n_traj <- 10 # 10000

  # Number of projections
  n_proj <- 10 # 1000
  
  full_data <- list("outbreak_data" = outbreak_data, "length_obs" = length_obs, "serial_interval" = serial_interval, 
                    "mu" = mu, "sigma" = sigma, "delta" = delta, "no_chunks" = no_chunks, 
                    "n_traj" = n_traj, "n_proj" = n_proj)
  
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

R_estimate <- function(observed_incidence, mu, sigma) {
  # Estimate R
  library(earlyR)
  # The R estimation
  R <- get_R(observed_incidence, si_mean = mu, si_sd = sigma, 
             max_R = 10)
  return(R)
}

projection <- function(full_data, obs_incidence, cutoff_time, proj_start, proj_end) {
  # Forecast incidence
  library(projections)

  # number of days that I hid data for i.e. number of days I'm projecting for
  hidden_days <- proj_end - cutoff_time
  
  # Call the R_estimate function
  R <- R_estimate(obs_incidence[1:cutoff_time, ], full_data$mu, full_data$sigma)
  
  # Project
  proj <- project(obs_incidence[1:cutoff_time, ], R = sample_R(R, 1000), si = full_data$serial_interval, 
                  n_sim = full_data$n_traj, n_days = hidden_days, R_fix_within = TRUE)
  return(proj)
}

# Do many projections
# multi_projection <- function(data) {
#   end <- length(data$daily_incidence) # when does the data end
#   hidden_days <- end - data$length_obs
#   
#   multi_projection <- array(NA, c(hidden_days, data$n_traj, data$n_proj)) # create array for storing each projection
#   
#   for (i in 1:data$n_proj){
#     single_projection <- projection(data$daily_incidence, data$length_obs, data$serial_interval, 
#                                     data$mu, data$sigma, data$n_traj)
#     multi_projection[, , i] <- single_projection   
#   }
#   return(multi_projection)
# }

hash_function <- function(projected_incidence) {
  library("digest")
  hash <- digest(projected_incidence)
  return(hash)
}

# reformat_dataset <- function(multi_projection, data) {
#   dataset <- data.frame(matrix(NA, nrow = 0, ncol = 3))
#   for (i in 1:data$n_proj){
#     for (j in 1:data$n_traj){
#       hash <- array(hash_function(multi_projection[ , j, i]), c(length(multi_projection[ , j, i]), 1))
#       time_seq <- seq(1, length(multi_projection[ , j, i]), 1)
#       timepoint <- array(time_seq, c(length(multi_projection[ , j, i]), 1))
#       hash_time <- cbind(hash, timepoint)
#       trajectory_data <- cbind(hash_time, multi_projection[ , j, i])
#       dataset <- rbind(dataset, trajectory_data)
#     }
#   }
#   colnames(dataset) <- c("hash", "days_since_data", "num_cases")
#   return(dataset)
# }

###########################
## The projection output ##
###########################

output <- function() {
  # Set seed for reproducibility
  set.seed(1)
  full_data <- full_data() # has all the original data I need 
  obs_incidence <- observed_incidence(full_data$outbreak_data) # has the incidence object for the outbreak
  combo_list <- split_data(full_data) # all the combinations of deltas that I want to project for
  
  for (i in 1:nrow(combo_list)) {
    cutoff_time <- full_data$delta * combo_list[i, 1] # the time at which observed data stops
    proj_end <- full_data$delta * combo_list[i, 2] # the time at which projection ends
    proj_start <- cutoff_time + 1 # the time at which projection starts
    
    cutoff_time <- full_data$delta * combo_list[1, 1] # the time at which observed data stops
    proj_end <- full_data$delta * combo_list[1, 2] # the time at which projection ends
    proj_start <- cutoff_time + 1 # the time at which projection starts
    
    proj_window <- projection(full_data, obs_incidence, cutoff_time, proj_start, proj_end)
  }
  
  # test <- multi_projection(full_data)
  # reformat_test <- reformat_dataset(test, full_data)
  # write.csv(test, file = "test.csv")
  return(proj_window)
}

