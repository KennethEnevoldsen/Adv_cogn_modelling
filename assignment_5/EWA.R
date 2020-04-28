#'@title Experience weighted attraction model
#'@description
#'
#'@param phi memory of old trials
#'@param delta weighting of forgone vs recieved
#'@param rho discounting of old trials
#'@param lambda consistency of choice with attraction (inverse heat)
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#' The actions of the agent
#'
#'@references
#'
#'@export
EWA <- function(n_agents=100, n_trials=20, n_tokens=20, pi, delta, rho, phi, lambda){
  
  # setup
  c <- array(0, c(n_agents, n_trials))
  n <- array(0, c(n_agents, n_trials))
  A <- array(0, c(n_agents, n_trials, n_tokens))
  p <- array(0, c(n_agents, n_trials, n_tokens))

  n[, 1] <- 1
  c[,1] <- 20
  A[,1,] <- 1
  
  for (t in 2:n_trials){
    for (a in 1:n_agents){
      n[a, t] <- (rho[a]*n[a, t-1]) + 1  # eq. 2.1 (experience updating)
      
      for (tok in 1:n_tokens){
        A[a, t, tok] <- (
          (phi[a]*n[a, t-1]*A[a, t-1, tok]) +  # prior attraction
            (delta[a] + ((1-delta[a])*(c[a, t-1] == tok))) *  # indicates whether the token was chosen
            ((((tok + sum(c[-a, t-1]))*pi)/n_agents)-tok)  # payoff for each possible contrib.
        )/ n[a, t]  # experience weighting
      }
      p[a, t, ] <- exp(lambda[a]*A[a, t,])/sum(exp(lambda[a]*A[a, t, ]))  # softmax
      c[a, t] <- extraDistr::rcat(1, p[a, t, ])
    }
  }
  internal_states=list(delta=delta, rho=rho, phi=phi, lambda=lambda)
  res <- list(choice=c, n=n, prob=p, pi=pi, internal_states=internal_states)
  return(res)
}


# test:
# EWA(n_agents=100, n_trials=20, n_tokens=20, pi=1.5, delta=runif(100, 0.2, 0.4), rho=runif(100, 0.2, 0.4), phi=runif(100, 0.2, 0.4), lambda=runif(100, 0.8, 2))
# n_agents=100; n_trials=20; n_tokens=20; pi=1.5; delta=runif(100, 0.2, 0.4); rho=runif(100, 0.2, 0.4); phi=runif(100, 0.2, 0.4); lambda=runif(100, 0.8, 2)
