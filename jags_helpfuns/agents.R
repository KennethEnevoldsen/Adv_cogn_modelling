
#'@title RB
#'@description
#'A random biased agent
#'
#'@param 
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'the reward and choice of the RB agent
#'
#'@references
#'
#'@export
RB <- function(payoff, theta){
  ntrials = nrow(payoff)
  b <- c(theta, 1-theta)
  choices <- array(0, c(ntrials))
  reward <- array(0, c(ntrials))
  
  for (t in 1:ntrials){
    # agent that chooses randomly between options 1 and 2 bias (to 1) = theta
    choices[t] <- extraDistr::rcat(1, b) # categorical distribution
    
    #  what reward does the agent get
    reward[t] <- payoff[t, choices[t]]
  }
  
  return(list(reward = reward, choices = choices))
}
# # Test
# biased_agent(ntrials = 3, theta = 0.5, payoff = cbind(c(0,0,0), c(1,1,1)))



#'@title Rescorla Wagner learning agent
#'@description
#'trivial from title
#'
#'@param payoff A list of payoffs
#'@param alpha The learning rate which should be between 0 and 1
#'@param beta The behavioural temperature, higher is more consistent. Should be below 0.
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'a list containing choice, reward, Q and p
#'
#'@references
#'(Rescorla et al., 1972)
#'
#'@export
rw <- function(payoff, alpha, beta){
  n_trials = nrow(payoff)
  choice <- array(0, c(n_trials))
  r <- array(0, c(n_trials))
  Q <- array(0, c(n_trials, 2))
  
  p <- array(0, c(n_trials, 2))
  
  # trial
  Q[1,1] <- 1
  Q[1,2] <- 1
  
  choice[1] <- extraDistr::rcat(1, c(0.5, 0.5))
  r[1] <- payoff[1, choice[1]]
  
  for (t in 2:n_trials){
    
    # learn
    Q[t,] <- Q[t-1,] + alpha*(r[t-1] - Q[t-1, ])
    Q[t, 2-(choice[t-1]-1)] <- Q[t-1, 2-(choice[t-1]-1)] # but only on the given choice
    
    # make choice
    p[t,] <-  exp(beta*Q[t,])/sum( exp(beta*Q[t,]))
    choice[t] <- rcat(1, p[t,])
    # get reward
    r[t] <- payoff[t, choice[t]]
  }
  return(list(choice = choice, reward = r, Q = Q, p = p))
}


