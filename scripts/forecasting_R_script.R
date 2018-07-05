# R script for forecasting an epidemic's trajectory

#########################################
## Functions for input into prediction ##
#########################################

# A discretised serial interval
si_disc <- function(mu, cv) {
  if(!require(epitrix)) stop("epitrix is missing")
  if(!require(distcrete)) stop("distcrete is missing")
  sim_si <- gamma_mucv2shapescale(mu, cv)
  si <- distcrete("gamma", shape = sim_si$shape, scale = sim_si$shape, w = 0, interval = 1)
  return(si)
}

# A simulated outbreak
sim_outbreak <- function(si, obs_time, R0) {
  if(!require(outbreaker)) stop("outbreaker is missing")
  # Simulate outbreak
  # set.seed(4)
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
  n_traj <- 10000 # 10000

  full_data <- list("outbreak_data" = outbreak_data, "serial_interval" = serial_interval, 
                    "mu" = mu, "sigma" = sigma, "delta" = delta, "no_chunks" = no_chunks, 
                    "n_traj" = n_traj)
  
  return(full_data)
}

# Listing the calibration and projection combinations
split_data <- function(outbreak_data) {
  # project for maximum 4 deltas
  # latest projection can be last chunk of no_chunks
  outbreak_data <- full_data()
  all_combos <- expand.grid(1:(outbreak_data$no_chunks - 1), outbreak_data$no_chunks)
  all_combos <- setNames(all_combos, c("calibration","projection"))
  # Want to keep combinations where calibration + projection <= no_chunks OR diff(projection, calibration) <= 4
  for (i in 1:nrow(all_combos)) {
    if ((all_combos$projection[i] - all_combos$calibration[i]) > 4) {
      all_combos$projection[i] <- all_combos$calibration[i] + 4
    }
  }
  return(all_combos)
}

##############################
## Functions for projection ##
##############################

R_estimate <- function(obs_incidence, serial_interval) {
  # Estimate R
  if(!require(earlyR)) stop("earlyR is missing")
  # The R estimation
  R <- get_R(obs_incidence, si = serial_interval, 
             max_R = 10)
  return(R)
}

projection <- function(full_data, obs_incidence, cutoff_time, proj_start, proj_end) {
  # Forecast incidence
  if(!require(projections)) stop("projections is missing")

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
# does the data look like it comes from the predictive probability distribution?
reliability <- function(x_t, pred) {
  if(!require(incidence)) stop("incidence is missing")
    reliability <- array(NA, dim = c(nrow(x_t))) # nrow(data) is how many days we have hidden data for
    for (i in 1:nrow(x_t)) {
      # lower_p <- (ppois((x_t - 1), mean(pred[i, ])) + 1e-7) # don't want to include the p into the boundary
      # upper_p <- (ppois((x_t + 1), mean(pred[i, ])) - 1e-7)
      # reliability[i] <- runif(1, min = lower_p, max = upper_p)
      
      reliability[i] <- ppois(x_t[i], mean(pred[i, ])) # cumulative_poisson(x_t[i], pred[i, ])
      # ppois calculates the cumulative probability distribution for time t: F_t
      # want the probability that there is x_t number of cases
    }
    # Next you apply a test of uniformity to these CDFs - chi-squared?
    # If they're uniform, the predictive distribution is the true data distribution
    
    return(reliability)
}

uniformity <- function(reliability) {
  if (sum(reliability) == 0) {
    uniformity <- NA
  } else{
    uniformity <- (chisq.test(reliability)$p.value)
  }
  return(uniformity)
}

ad_test <- function(reliability) {
  # devtools::install_github("bbolker/ADmarsaglia")
  if(!require(ADmarsaglia)) stop("ADmarsaglia is missing")
  ordered_array <- reliability[order(reliability)] # already a CDF - but need to order smallest to largest
  ad_score <- ADtest(ordered_array)
  # ad_pvalue <- AD_pval(length(reliability), ad_score)
  # shapiro <- array(NA, dim = c(length(reliability))) # empty array
  # 
  # for (i in 1:length(reliability)) {
  #   shapiro[i] <- (((2 * i) - 1) / length(reliability)) * (log(punif(ordered_array[i])) + 
  #                 log(punif(1 - ordered_array[length(reliability) + 1 - i])))
  # }
  # 
  # sum_shapiro <- sum(shapiro)
  # ad_score <- -length(reliability) - sum_shapiro
  # if (ad_score == Inf) {
  #   ad_score <- 1000
  # }
  return(ad_score$p.value)
}


## Sharpness

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
    # don't want to divide by 0 but sharpness is relative - move the projections a bit
    pred_ind <- as.vector(pred[i, ] + 1)

    sharpness[i] <- 1 - (MADM(pred_ind) / median(pred_ind)) # median is the forecast median
  }
  return(sharpness)
}

# sum_sharpness <- function(pred){
#   # array for storing daily sharpness
#   # sum up the number of cases in a prediction window for each projection and avoid it being 0
#   # don't want to divide by 0 but sharpness is relative - move the projections a bit
#   sum_pred <- array(NA, dim = c(ncol(pred)))
#   for (i in 1:ncol(pred)){
#     sum_pred[i] <- as.vector(sum(pred[ , i]) + 1)
#   }
#   sharpness <- 1 - (MADM(sum_pred) / median(sum_pred))
#   return(sharpness)
# }



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

# Aggregated RMSE for a given time window
# Feed in the $n of the incidence object
# sum_rmse <- function(x_t, pred){
#   sum_pred <- array(NA, dim = c(ncol(pred)))
#   for (i in 1:ncol(pred)){
#     sum_pred[i] <- as.vector(sum(pred[ , i]))
#   }
#   rmse <- sqrt(mean_squared_diff(x_t, sum_pred))
#   return(rmse)
# }



#################################################
## The projection and prediction metric output ##
#################################################

output <- function() {
  if(!require(dplyr)) stop("dplyr is missing")
  if(!require(incidence)) stop("incidence is missing")
  if(!require(digest)) stop("digest is missing")
  
  # Set seed for reproducibility
  # set.seed(1)
  sim_gen_data <- full_data() # has all the original data I need
  
  # Create a directory into which I store this simulation
  sim_hash <- digest(sim_gen_data) # unique identifier for simulation
  dir.create(paste("/home/evelina/Development/forecasting/simulations/ebola_test_4/", sim_hash, sep = ""))
  
  # Set this directory as the working directory so files end up in the right place
  setwd(paste("/home/evelina/Development/forecasting/simulations/ebola_test_4/", sim_hash, sep = ""))
  
  # Projections and prediction metrics
  # Incidence object for the outbreak
  linelist <- sim_gen_data$outbreak_data
  obs_incidence <- incidence(linelist$onset, interval = 1, last_date = (sim_gen_data$delta * sim_gen_data$no_chunks))
  # All the combinations of deltas that I want to project for
  combo_list <- split_data(sim_gen_data) 
  
  # Save plot of incidence curve
  pdf("incidence.pdf", width = 7, height = 5)
    print(plot(obs_incidence))
  dev.off()    
  
  # A for loop for doing all the projections and calculating prediction metrics
  for (i in 1:nrow(combo_list)) {
    cutoff_time <- sim_gen_data$delta * combo_list[i, 1] # the time at which observed data stops
    window_n <- combo_list[i, 2] - combo_list[i, 1] # number of projection windows
    proj_end <- sim_gen_data$delta * combo_list[i, 2] # the time at which projection ends
    proj_start <- cutoff_time + 1 # the time at which projection starts
    

    # if (i == 7) {
    #   print("We should have NAs here")
    # }
    
    proj_window <- projection(sim_gen_data, obs_incidence, cutoff_time, proj_start, proj_end) # projection for the time window
    
    pdf(paste("projection", i, ".pdf", sep = ""), width = 7, height = 5)
      print(plot(obs_incidence) %>% add_projections(proj_window))
    dev.off()
    
    proj_rel <- reliability(obs_incidence[proj_start:proj_end, ]$counts, proj_window) # daily reliability
    proj_bias <- bias(obs_incidence[proj_start:proj_end, ]$counts, proj_window) # daily bias
    proj_sharp <- sharpness(proj_window)
    proj_rmse <- rmse(obs_incidence[proj_start:proj_end, ]$counts, proj_window)
    
    # Start up arrays for storing the metrics for all projection windows
    proj_rel_test <- array(NA, dim = c(window_n)) # reliability
    proj_bias_mean <- array(NA, dim = c(window_n)) # bias
    proj_sharp_mean <- array(NA, dim = c(window_n)) # sharpness
    proj_rmse_mean <- array(NA, dim = c(window_n)) # root-mean-square-error
    proj_group <- array(NA, dim = c(window_n)) # identifies what types of projections the window contains
    
    # Calculate prediction metrics by projection window    
    for (j in 1:window_n){
      window_end <- sim_gen_data$delta * j # the time at which projection window ends
      window_start <- window_end - (sim_gen_data$delta - 1) # the time at which projection window starts
      
      # proj_rel_test[j] <- uniformity(proj_rel[window_start:window_end]) # window-specific reliability
      # proj_bias_score[j] <- mean(proj_bias[window_start:window_end])
      # proj_sharp[j] <- sum_sharpness(proj_window[window_start:window_end, ])
      # proj_rmse_sum[j] <- rmse(obs_incidence[window_start:window_end, ]$n, proj_window[window_start:window_end])
      
      proj_rel_test[j] <- ad_test(proj_rel[window_start:window_end]) # window-specific reliability
      proj_sharp_mean[j] <- mean(proj_sharp[window_start:window_end]) # window-specific sharpness
      proj_bias_mean[j] <- mean(proj_bias[window_start:window_end]) # window-specific bias
      proj_rmse_mean[j] <- mean(proj_rmse[window_start:window_end]) # window-specific RMSE

      if (sum(proj_window[window_start:window_end, ]) == 0) {
        proj_group[j] <- 1 # the projections only contain zeroes
      } else if (any(colSums(proj_window[window_start:window_end, ]) == 0)) {
        proj_group[j] <- 2 # the projections contain some zeroes
      } else {
        proj_group[j] <- 3 # there are no projections that are just zeroes
      }
    }

    # Window-specific metrics
    # proj_rel_test <- array(NA, dim = c(combo_list[i, 2] - combo_list[i, 1]))
    # proj_rmse_sum <- array(NA, dim = c(combo_list[i, 2] - combo_list[i, 1]))
    # proj_bias_score <- array(NA, dim = c(combo_list[i, 2] - combo_list[i, 1]))
    # for (j in 1:(combo_list[i, 2] - combo_list[i, 1])){
    #   window_end <- sim_gen_data$delta * j
    #   window_start <- (window_end - sim_gen_data$delta) + 1
    #   proj_rel_test[j] <- uniformity(proj_rel[window_start:window_end])
    #   proj_rmse_sum[j] <- sum(proj_rmse[window_start:window_end], na.rm = TRUE)
    #   proj_bias_score[j] <- mean(proj_bias[window_start:window_end])
    # }
    
    # make an array for storing the prediction metrics for a given projection
    proj_metrics <- data.frame(dataset = sim_hash,
                               cali_window_size = cutoff_time,
                               no_cali_cases = obs_incidence[1:cutoff_time, ]$n,
                               proj_window_size = sim_gen_data$delta,
                               proj_window_no = c(1:(combo_list[i, 2] - combo_list[i, 1])),
                               disease = "ebola",
                               reliability = proj_rel_test, 
                               sharpness = proj_sharp_mean,
                               bias = proj_bias_mean,
                               rmse = proj_rmse_mean,
                               pred_type = proj_group)

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
    output()
    print(paste("Finished simulation", i, sep = " "))
  }
  return(print("Finished all projections"))
}

# multi_output(1)