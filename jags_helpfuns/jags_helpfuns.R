
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
bandit_payoff(100, probs = c(.3, .8), reward = c(1, 1.5))


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
sim_fit_compare <- function(payoff_gen_fun, agent_funs, params_to_save, model_filepath, n_sim = 1){
  if (n_sim > 1){
    for (i in 1:n_sim){
      print(paste("Currenty running simulation: ", i, " out of ", n_sim, sep = ""))
      tmp = sim_fit_compare(payoff_gen_fun, agent_funs, params_to_save, model_filepath, n_sim = 1)
      tmp$n_sim = i
      if (i == 1){
        res <- tmp
      } else {
        res <- rbind(res, tmp)
      }
    }
    return(res)
  }
  res = array(0, c(length(agent_funs), length(agent_funs)))
  res = as.data.frame(res)
  colnames(res) = names(agent_funs)
  res$model_fitted_to_data = rep("NA", nrow(res))
  
  payoff <- eval(parse(text = payoff_gen_fun))
  
  for (a in 1:length(agent_funs)){
    nam = names(agent_funs[a])
    print(paste("Currently similating and fitting to agent: ", nam, sep = ""))
    sim_dat <- eval(parse(text = agent_funs[[a]])) # simulate data
    
    # fit each model to data
    for (i in 1:length(params_to_save)){
      samples <- jags.parallel(data = list(choice = sim_dat$choice, n_trials = length(sim_dat$choice), r = sim_dat$reward), 
                               inits = NULL, 
                               parameters.to.save = params_to_save[[i]], 
                               model.file = model_filepath[[i]],
                               n.chains = 4, n.iter = 3000, n.burnin = 1000)
      res[i, nam] <-  samples$BUGSoutput$DIC
      res$model_fitted_to_data[i] <- names(params_to_save[i])
    }
  }
  return(res)
}