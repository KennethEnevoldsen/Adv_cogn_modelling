
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


#'@title Bandit payoff
#'@description
#'Generate a payoff for an n bandit task, where n is the number of bandits
#'n is derived from the length of probs and reward
#'
#'@param props the probability of reward, can be given as a list should be the same length as reward
#'@param reward the reward for a given prob
#'@param loss the expected loss if NULL (default) assumed to be 0
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'a matrix of size n_trials x n
#'
#'@references
#'
#'
#'@export
bandit_payoff <- function(n_trials = 100, probs = c(0.3, 0.8), reward = c(1, 1.5), loss = NULL){
  if (length(probs) != length(reward)){
    stop("length(probs) != length(reward), which the function assumes, make sure it is the case")
  } 
  if (is.null(loss)){
    loss <- rep(0, length(probs))
  } else if (length(loss) != length(probs)){
    stop("length(probs) != length(loss), which the function assumes, make sure it is the case")
  }
  
  res =  array(0, c(n_trials, length(probs)))
  for (i in 1:length(probs)){
    res[,i] <- rbinom(n_trials, 1, probs[i])*reward[i]
    res[,i][res[,i] == 0] = loss [[i]]
  }
  return(res)
}
# bandit_payoff(100, probs = c(.3, .8), reward = c(1, 1.5))


#'@title simulate fit and compare your jags models and agents
#'@description
#'
#'
#'@param todo
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'a dataframe
#'
#'@references
#'
#'
#'@export
simulate_fit <- function(gen_fun, data_to_fit, model_filepath, params_to_save, save_samples = F, n_sim = 1){
  
  if (n_sim > 1){
    for (i in 1:n_sim){
      print(paste("Currenty running simulation: ", i, " out of ", n_sim, sep = ""))
      tmp = simulate_fit(gen_fun, data_to_fit, model_filepath, params_to_save, save_samples,  n_sim = 1)
      tmp$n_sim = i
      if (i == 1){
        res <- tmp
      } else {
        res <- rbind(res, tmp)
      }
    }
    return(res)
  }

  if (length(unique(names(gen_fun))) != length(gen_fun)){
    stop("The entries in the gen_fun list must be named with unique names")}
  if (length(unique(names(model_filepath))) != length(model_filepath)){
    stop("The entries in the model_filepath list must be named with unique names")}
  
  res_l = NULL
  for (a in 1:length(gen_fun)){
    nam <- names(gen_fun[a])
    print(paste("Currently similating and fitting to the generative model: ", nam, sep = ""))
    sim_dat <- eval(parse(text = gen_fun[[a]])) # simulate data
    
    # make dataframe to save in
    res <- dplyr::tibble("model_generating_the_data" = rep(nam, length(params_to_save)), 
                  "DIC" = NA,
                  "model_fitted_to_data" = NA)
    if (isTRUE(save_samples)){
      res$samples <- NA
      res$true_params <- rep(list(sim_dat$start_params), length(params_to_save))
    }
    
    # fit each model to data
    for (i in 1:length(params_to_save)){
      # cat("\nlength params to save: ", length(params_to_save))
      samples <- R2jags::jags.parallel(data = eval(parse(text = data_to_fit[[i]])), 
                               inits = NULL, 
                               parameters.to.save = params_to_save[[i]], 
                               model.file = model_filepath[[i]],
                               n.chains = 4, n.iter = 5000, n.burnin = 1000, n.thin = 2)
      res$DIC[i] <-  samples$BUGSoutput$DIC
      res$model_fitted_to_data[i] <- names(model_filepath[i])
      if (isTRUE(save_samples)){
        res$samples[i] <- list(samples)
      }
    }
    res_l[[a]] <- res
    rm(sim_dat)
  }
  
  res <- res_l %>% 
    do.call(rbind, .) 
  return(res)
}

               
