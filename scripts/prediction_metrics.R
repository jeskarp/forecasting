# Metrics for evaluating projection performance
# Janetta Skarp


## Average residual
# x_t is the true daily incidence data
# pred is the projected daily incidence data

residual_difference <- function(x_t, pred){
  # Calculate difference between data point and sample
  difference <- array(NA, dim = c(length(pred)))
  for (i in 1:length(pred)){
    difference[i] <- x_t - pred[i]
  }
  return(difference)
}

residual <- function(x_t, pred) {
  residual <- array(NA, dim = c(length(x_t)))
  for (i in 1:length(x_t)) {
    residual[i] <- mean(residual_difference(x_t[i], as.vector(pred[i, ])))
  }
  return(residual)
}



## Mean-squared error
# x_t is the true daily incidence data
# pred is the projected daily incidence data

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

# Calculates MSE for each forecasted timepoint
# x_t is the true incidence data, pred is the projected incidence data
mse <- function(x_t, pred){
  mse <- array(NA, dim = c(length(x_t)))
  for (i in 1:length(x_t)){
    mse[i] <- (mean_squared_diff(x_t[i], as.vector(pred[i, ]))) / (x_t[i] + 1) 
    # add +1 to MSE because don't want to divide by 0
  }
  return(mse)
}



## Sharpness
# pred is the projected daily incidence data

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



## Bias
# x_t is the true daily incidence data
# pred is the projected incidence data

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