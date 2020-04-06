
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
  
  return(list(reward = reward, 
              choices = choices, 
              start_params = list(theta = theta)))
}
# # Test
# biased_agent(ntrials = 3, theta = 0.5, payoff = cbind(c(0,0,0), c(1,1,1)))



#'@title Rescorla Wagner learning agent
#'@description
#'trivial from title
#'
#'@param payoff A list of payoffs
#'@param alpha The learning rate which should be between 0 and 1
#'@param beta The behavioural temperature, higher is more consistent. Should not be below 0.
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
  Q_update <- array(0, c(n_trials, 2))
  p <- array(0, c(n_trials, 2))
  
  # trial
  Q[1,1] <- 1
  Q[1,2] <- 1
  
  choice[1] <- extraDistr::rcat(1, c(0.5, 0.5))
  r[1] <- payoff[1, choice[1]]
  
  for (t in 2:n_trials){
      for (k in 1:2){
        # update utility Q for chosen option with reward on last trials
        Q_update[t, k] <- Q[t-1, k] + alpha*(r[t-1] - Q[t-1, k])
        Q[t, k] <- ifelse(k==choice[t-1], Q_update[t, k], Q[t-1, k])
        exp_p[t, k] <- exp(beta*Q[t,k])
      }
      for (k in 1:2){
        p[t,k] <- exp_p[t,k]/sum(exp_p[t,])
      }
    # make choice
    p[t,] <-  exp(beta*Q[t,])/sum( exp(beta*Q[t,]))
    choice[t] <- rcat(1, p[t,])
    # get reward
    r[t] <- payoff[t, choice[t]]
  }
  return(list(choice = choice, 
              reward = r, 
              Q = Q, 
              p = p,
         start_params = list(alpha = alpha, beta = beta)))
}

#'@title Rescorla Wagner learning agent
#'@description
#'trivial from title
#'
#'@param payoff A list of payoffs
#'@param alpha The learning rate which should be between 0 and 1
#'@param beta The behavioural temperature, higher is more consistent. Should not be below 0.
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'a list containing choice, reward, choice kernel and probability of choice, and starting parameters
#'
#'@references
#'
#'
#'@export
choice_kernel <- function(payoff, alpha, beta){
  n_trials = nrow(payoff)
  choice <- array(0, c(n_trials))
  r <- array(0, c(n_trials))
  ck <- array(0, c(n_trials, 2)) # choice kernel
  ck_chosen <- array(0, c(n_trials, 2)) 
  ck_unchosen <- array(0, c(n_trials, 2)) 
  p <- array(0, c(n_trials, 2))
  exp_p <- array(0, c(n_trials, 2))
  
  # trial1
  ck[1,1] <- 1
  ck[1,2] <- 1
  
  choice[1] <- extraDistr::rcat(1, c(0.5, 0.5))
  r[1] <- payoff[1, choice[1]]
  
  
  for (t in 2:n_trials){
    for (k in 1:2){
      # learn
      ck_chosen[t,k] <- ck[t-1,k] + alpha*(1-ck[t-1,k])
      ck_unchosen[t,k] <- ck[t-1,k] + alpha*(0-ck[t-1,k])
      
      ck[t, k] <- ifelse(k == choice[t-1], ck_chosen[t, k], ck_unchosen[t, k])
      
      exp_p[t, k] <- exp(beta * ck[t, k])
    }
    for (k in 1:2){
      p[t, k] <- exp_p[t, k] / sum(exp_p[t,])
    }
    
    choice[t] <- extraDistr::rcat(1, p[t,])
    # get reward
    r[t] <- payoff[t, choice[t]]
  }
  
  return(list(choice = choice, 
              reward = r, 
              ck = ck, 
              p = p,
              start_params = list(alpha = alpha, beta = beta)))
}




#'@title PVL delta Agent
#'@description
#'
#'
#'@param a The learning rate of the agent 0 <= a <= 1
#'@param A The shape parameter  0 <= A <= 1
#'@param w The loss aversion parameter. 0 <= w <= 5
#'@param c response consistency. 0 <= c <= 5
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'
#'
#'@references
#'
#'@export
pvl_delta <- function(payoff, a, A, w, theta){
  if (is.list(payoff)){
    payoff = payoff$payoff
  }
  
  
  n_trials <- nrow(payoff) # number of trials
  n_decks <- ncol(payoff) # number of decks
  
  choice <- rep(0, n_trials)
  r <- rep(0, n_trials)
  u <- array(0, c(n_trials, n_decks))
  ev <- array(0, c(n_trials, n_decks))
  p <- array(0, c(n_trials, n_decks))
  
  tmp_p = array(0, c(n_trials, n_decks))
  
  p[1,] <- rep(1/n_decks, n_decks)
  choice[1] <- extraDistr::rcat(1, p[1,])
  r[1] <- payoff[1, choice[1]]
  
  for (t in 2:n_trials){ # for each trial
    for (d in 1:n_decks){ # for each deck
      
      # calculate subjective utility
      u[t,d] <- ifelse(r[t-1] >= 0, 
                       abs(r[t-1])^A,
                       -w * abs(r[t-1])^A
      )
      
      ev[t,d] <- ifelse(choice[t-1] == d, 
                        ev[t-1, d] + a * (u[t, d] - ev[t-1, d]),
                        ev[t-1, d]
      )
      # for softmax
      tmp_p[t, d] <- exp(theta * ev[t, d])
      if (is.infinite(tmp_p[t, d])){
        if (tmp_p[t, d] > 0){
          tmp_p[t, d] = exp(500)
        } 
      } 
    }
    
    # update prop
    for (d in 1:n_decks){
      p[t,d] <- tmp_p[t,d]/sum(tmp_p[t,1:n_decks])
    }
    choice[t] <- extraDistr::rcat(1, p[t,]) # categorical distribution
    
    #  what reward does the agent get
    r[t] <- payoff[t, choice[t]]
  }
  
  return(list(reward = r, 
              choices = choice, 
              ev = ev,
              p = p,
              u = u,
              payoff = payoff,
              start_params = list(a = a, A = A, w = w, theta = theta)))
}


#'@title VSE Agent
#'@description
#'
#'
#'@param reward a vector containing reward or a list containg both loss and reward
#'@param loss  a vector containing loss or NULL if reward contains both
#'@param alpha learning rate (bound between 0 and 1) - 0 is no learning
#'@param beta inverse heat. degree of exploration - beta = 0 is completely random exploration
#'@param delta decay parameter (bounded between 0 and 1) - 0 mean higher reliance on recent outcomes
#'@param theta value sensitivity, similar to the PVL delta model (bounded between 0 and 1)
#'@param phi exploration bonus (unbounded)
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'a list containing the choice and internal states of the agent
#'
#'@references
#'Ligneul, 2019
#'@export
vse <- function(reward, loss = NULL, alpha, theta, delta, phi, beta){
  if (is.list(reward)){
    loss = reward$loss
    reward = reward$reward
  } 
  if (is.null(loss) & !is.list(reward)) {
    stop("loss is NULL but reward is not a list. Therefore the model lacks a reward")
  }
  n_trials <- nrow(reward) # number of trials
  n_decks <- ncol(reward) # number of decks
  
  choice <- rep(0, n_trials)
  v <- array(0, c(n_trials))
  exploit <- array(0, c(n_trials, n_decks))
  explore <- array(0, c(n_trials, n_decks))
  p <- array(0, c(n_trials, n_decks))
  tmp_p = array(0, c(n_trials, n_decks))
  
  p[1,] <- rep(1/n_decks, n_decks)
  choice[1] <- extraDistr::rcat(1, p[1,])
  
  for (t in 2:n_trials){ # for each trial
    
    # calculate 'utiity'
    v[t] <- reward[t, choice[t-1]]^theta - loss[t-1, choice[t-1]]^theta
    
    for (d in 1:n_decks){ # for each deck
      exploit[t, d] <- ifelse(choice[t-1] == d,
                              exploit[t-1, d] * delta + v[t],
                              exploit[t-1, d] * delta)
      explore[t, d] <- explore[t-1, d] + alpha * (phi - explore[t-1, d])
      
      # for softmax
      tmp_p[t, d] <- exp(beta * (exploit[t, d] + explore[t, d]))
      if (is.infinite(tmp_p[t, d])){
        if (tmp_p[t, d] > 0){
          tmp_p[t, d] = exp(500)
        } 
      } 
    }
    
    # update prop
    for (d in 1:n_decks){
      p[t,d] <- tmp_p[t,d]/sum(tmp_p[t,1:n_decks])
    }
    choice[t] <- extraDistr::rcat(1, p[t,]) # categorical distribution
    
  }
  
  return(list(reward = reward,
              loss = loss,
              choices = choice, 
              exploit = exploit,
              explore = explore,
              p = p, 
              v = v,
              start_params = list(alpha = alpha,
                                  theta = theta, 
                                  delta = delta, 
                                  phi = phi, 
                                  beta = beta)))
}

#'@title ORL Agent
#'@description
#'
#'
#'@param a_rew learning rate for rewards
#'@param a_pun = learning rate for punishments
#'@param w_f influnce of win/loss frequency on valence (relative to value)
#'@param w_p influence of perseverance on valence (rel. to value)
#'@param decay decay for perseverance
#'@param theta softmax heat
#'
#'@author
#'K. Enevoldsen
#'
#'@return 
#'
#'
#'@references
#'
#'@export
orl <- function(payoff, a_rew, a_pun, w_f, w_p, decay) {
  n_trials <- nrow(payoff)
  n_decks <- ncol(payoff)
  c <- n_decks - 1
  
  
  choice <- array(0,c(n_trials))
  r <- array(0,c(n_trials))
  # expected value
  ev <- array(0,c(n_trials, n_decks))
  # expected value update
  ev_update <- array(0,c(n_trials, n_decks))
  # expected frequency of reward
  ef <- array(0, c(n_trials, n_decks))
  # (temporary use) updated expected frequency
  ef_chosen <- array(0, c(n_trials, n_decks))
  # (temp) updated expected frequency for not chosen decks
  ef_not_chosen <- array(0, c(n_trials, n_decks))
  sign_x <- array(0, c(n_trials))
  # perseverance
  pers <- array(0, c(n_trials, n_decks))
  # valence
  v <- array(0, c(n_trials, n_decks))
  # probability for choices
  p <- array(0,c(n_trials, n_decks))
  tmp_p = array(0, c(n_trials, n_decks))
  
  
  # choice first round
  p[1,] <-  c(0.25, 0.25, 0.25, 0.25)
  choice[1] <- extraDistr::rcat(1, p[1,])
  r[1] <- payoff[1, choice[1]]
  
  for (t in 2:n_trials){ # for each trial
    
    sign_x[t-1] <- ifelse(r[t-1] >= 0, 1, -1)
    if (r[t-1]==0){sign_x[t-1] <- 0}
    
    for (d in 1:n_decks){ # for each deck
      
      #estimated value
      ev_update[t, d] <- ifelse(r[t-1] >= 0, 
                                ev[t-1, d] + a_rew*(r[t-1] - ev[t-1, d]),
                                ev[t-1, d] + a_pun*(r[t-1] - ev[t-1, d])
      )            
      ev[t, d] <- ifelse(choice[t-1] == d, 
                         ev_update[t, d],
                         ev[t-1, d]
      )
      
      #estimated freq
      ef_chosen[t, d] <- ifelse(r[t-1] >= 0, 
                                ef[t-1, d] + a_rew*(sign_x[t-1] - ef[t-1, d]),
                                ef[t-1, d] + a_pun*(sign_x[t-1] - ef[t-1, d])
      )
      ef_not_chosen[t, d] <- ifelse(r[t-1] >= 0, 
                                    ef[t-1, d] + a_pun*(-(sign_x[t-1]/c) - ef[t-1, d]),
                                    ef[t-1, d] + a_rew*(-(sign_x[t-1]/c) - ef[t-1, d])
      )       
      ef[t,d] <- ifelse(choice[t-1] == d, 
                        ef_chosen[t, d],
                        ef_not_chosen[t, d]
      )
      
      # perseverence
      pers[t, d] <- ifelse(d == choice[t-1], 
                           1/(1+decay),
                           pers[t-1, d]/(1+decay)
      )
      
      # valence
      v[t, d] <- ev[t, d] + ef[t, d] * w_f + pers[t, d]*w_p
      
      # for softmax
      tmp_p[t, d] <- exp(v[t, d])
      if (is.infinite(tmp_p[t, d])){
        if (tmp_p[t, d] > 0){
          tmp_p[t, d] = exp(500)
        } 
      } 
      
    }
    
    # update prop
    for (d in 1:n_decks){
      p[t,d] <- tmp_p[t, d]/sum(tmp_p[t,1:n_decks])
    }
    
    
    choice[t] <- extraDistr::rcat(1, p[t,1:n_decks]) # categorical distribution
    
    #  what reward does the agent get
    r[t] <- payoff[t, choice[t]]
  }
  
  return(list(payoff = payoff,
              reward = r, 
              choices = choice, 
              ev = ev,
              p = p,
              ef = ef,
              pers = pers,
              v = v,
              start_params = list(a_rew = a_rew, a_pun = a_pun, w_f = w_f, w_p = w_p, decay = decay)))
}

