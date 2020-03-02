
#'@title get Maximum a posterior (MAP) from jags samples
#'@description
#'trivial from title
#'
#'@param param_samples Samples derived from jags
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'Maximum a posterior
#'
#'@references
#'
#'@export
jag_map <- function(param_samples){
  param_samples = as.numeric(param_samples)
  dens <- density(param_samples)
  return(dens$x[dens$y == max(dens$y)])
}

#'@title likelihood_with_x
#'@description
#'
#'
#'@param sample Samples derived from jags
#'@param true_state the value of the truth
#'@param epsilon
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'the probability to sample the truth +/- epsilon from the sample distribution
#'
#'@references
#'
#'@export
likelihood_with_x <- function(sample, true_state, epsilon = 0.01){
  min_t <- true_state - epsilon
  max_t <- true_state + epsilon
  # the probability to sample the truth +/- epsilon from the sample distribution
  return(length(sample[sample > min_t & sample < max_t])/length(sample)) 
}