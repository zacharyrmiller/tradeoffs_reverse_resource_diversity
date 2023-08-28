### zrmiller@illinois.edu | zachary.miller@yale.edu ###

### Functions for model simulation and visualization

##### Functions for C-R dynamics with and without stochasticity #####

full_CR_dynamics <- function(time, state, pars){
  
  ### deterministic (ODE) dynamics for full CR system (eq. 1-3 in Posfai et al.)
  
  # define model parameters
  A <- pars$A # matrix of resource utilization rates
  s <- pars$s # resource supply rates
  mu <- pars$mu # resource degradation rates
  v <- pars$v # resource values
  del <- pars$del # consumer death rates
  
  p <- pars$p # number of resources
  m <- pars$m # number of consumers
  
  r_function <- pars$r # functional response (shared for all resources)
  
  # define state variables
  cs <- state[1:p] # resource concentrations
  ns <- state[(p+1):(p+m)] # consumer abundances
  
  # compute rates of change
  dc <- s - diag(r_function(cs)) %*% (t(A) %*% ns) - mu * cs
  dn <- diag(ns) %*% (A %*% (v * r_function(cs)) - del)
  
  return(list(c(dc, dn)))
}


consumer_dynamics <- function(time, state, pars){
  
  ### simulate simplified ODE dynamics for consumers only (eq. 4 in Posfai et al.)
  
  # define model parameters
  A <- pars$A # matrix of resource utilization rates
  s <- pars$s # resource supply rates
  
  # define state variables
  n <- state # consumer abundances
  
  # compute rates of change
  dn <- diag(n) %*% (A %*% diag(1 / as.vector(t(A) %*% n)) %*% s - 1)
  
  return(list(as.vector(dn)))
}


stochastic_birth_death_dynamics <- function(state, 
                                pars, # a list of pars, or a flag for neutral dynamics
                                time = NA, # a vector of times at which births/deaths occur (return abundances for all times)
                                until = NA # a number of extinctions to simulate until (return time(s) to extinction)
){
  
  ### simulate stochastic birth-death process corresponding to consumer-only dynamics
  
  if(identical(pars, "neutral")){
    
    neutral <- TRUE
    
    g <- 1 # growth rates always equal for all consumers
  }else{
    
    neutral <- FALSE
    
    # define model parameters
    A <- pars$A # matrix of resource utilization rates
    s <- pars$s # resource supply rates
  }
  
  # define initial state variables
  n <- state # consumer abundances
  
  m <- length(n) # number of consumers
  S <- sum(n) # total number of individuals (= total resource inflow)
  
  # simulate one stochastic trajectory
  if(!identical(time, NA)){ # if time is supplied, simulate dynamics at all times
    
    out <- matrix(nrow = length(time), ncol = m) # initialize matrix of times x species
    
    for(i in 1:length(time)){
      
      this_step <- single_birth_death_step(n, m, A, s, S, neutral) # choose one individual to die and one to reproduce, according to specified dynamics
      
      n <- this_step$n # update consumer abundances
      
      # record new abundances
      out[i, ] <- n
    }
    
    # prepare output for return: add times to output and label columns
    out <- cbind(time, out)
    colnames(out) <- c("time", 1:m)
    
  }else{ # if time is not supplied, simulate until # of extinctions specified by "until" arg
    
    # initialize counters and output
    step_counter <- 0
    extinction_counter <- 0
    out <- c() # a vector of extinction times
    
    while(extinction_counter < until){
      
      step_counter <- step_counter + 1 # update time
      
      this_step <- single_birth_death_step(n, m, A, s, S, neutral) # choose one individual to die and one to reproduce, according to specified dynamics
      
      n <- this_step$n # update consumer abundances
        
      if(this_step$extinction){
        
        out <- append(out, step_counter) # add extinction to output
        extinction_counter <- extinction_counter + 1
      }
    }
  }
  
  return(out)
}


single_birth_death_step <- function(n, m, A, s, S, neutral){
  
  ### execute one step of stochastic birth-death dynamics, by choosing one individual to die and one to reproduce 
  ### dynamics can be neutral or reduced C-R, according to "neutral" arg
  
  extinction <- FALSE # flag indicating whether an extinction occured at this step
  
  # sample one individual to die (equally likely)
  which_dies <- sample(1:m, 1, prob = c(n / S))
  n[which_dies] <- n[which_dies] - 1
  
  # check if an extinction occured
  if(n[which_dies] == 0) extinction <- TRUE
  
  # compute growth rates
  if(!neutral) g <- A %*% (s / (t(A) %*% n))
  
  # sample one species to reproduce (proportional to growth rates)
  which_new <- sample(1:m, 1, prob = (g * n) / sum(g * n))
  n[which_new] <- n[which_new] + 1
  
  return(list(n = n, extinction = extinction))
}


stochastic_birth_death_immigration_dynamics <- function(state,
                                   pars,
                                   mutation_rate,
                                   max_steps,
                                   sample_every = 10^5){
  
  ### simulate stochastic birth-death-immigration process corresponding to consumer-only dynamics
  ### this function runs a set number of simulation steps, sampling abundances at specified intervals
  
  if(identical(pars, "neutral")){
    
    neutral <- TRUE
    
    g <- 1 # growth rates always equal
  }else{
    
    neutral <- FALSE
    
    # define model parameters
    A <- pars$A # matrix of resource utilization rates
    s <- pars$s # resource supply rates
  }
  
  # define initial state variables
  n <- state # consumer abundances
  
  m <- length(n) # number of consumers
  p <- length(s)
  S <- sum(n) # total number of individuals (= total resource inflow)
  
  out <- list()
  
  for(i in 1:max_steps){
    
    if((i %% sample_every) == 0) out <- append(out, list(n)) # record abundances at specified intervals
    
    # sample one individual to die (equally likely)
    which_dies <- sample(1:m, 1, prob = c(n / S))
    n[which_dies] <- n[which_dies] - 1
    
    if(n[which_dies] == 0){ # if a species goes extinct, update parameters accordingly
      
      m <- m - 1 # decrease number of consumers 
      if(!neutral) A <- A[-which_dies, , drop = FALSE] # remove parameters from A
      n <- n[-which_dies] # remove from vector of consumer abundances (to prevent this vector from growing too large)
    }
    
    if(runif(1) < mutation_rate){ # if a mutation/immigration event occurs, update parameters
      
      m <- m + 1 # increase number of consumers
      if(!neutral) A <- rbind(A, sample_from_simplex(p) + rnorm(p, 0, pars$sigma)) # sample a new species with random metabolic strategy
      n <- c(n, 1) # append new species (one individual) to vector fo consumber abundances
    }else{
      
      # compute growth rates
      if(!neutral) g <- A %*% (s / (t(A) %*% n))
      
      # sample one species to reproduce (proportional to growth rates)
      which_new <- sample(1:m, 1, prob = (g * n) / sum(g * n))
      n[which_new] <- n[which_new] + 1
    }
  }
  
  # organize output into a data frame
  out <- lapply(out, as.data.frame)
  out <- bind_rows(out, .id = "Name")
  out[] <- lapply(out, as.numeric)
  colnames(tmp) <- c("timepoint", "abuns")
  
  return(out)
}


##### Auxilliary functions #####

generate_tradeoff_matrix <- function(p, m){
  
  ### generate a m x p matrix where each row is a random uniform sample from the p - 1 simplex
  
  A <- matrix(rexp(p * m), nrow = m, ncol = p)
  A <- diag(1 / rowSums(A)) %*% A # rows sum to one
  
  return(A)
}


sample_from_simplex <- function(k){
  
  ### generate a random uniform sample from the k simplex
  
  out <- rexp(k)
  return(out / sum(out))
}


test_inhull <- function(A, s, simplex = TRUE){
  
  ### check if a point s is in the convex hull of the rows of A
  
  if(simplex){ # if s and A are already in the simplex, drop one coordinate
    
    A <- A[, -1, drop = FALSE]
    s <- s[-1]
  }
  
  s <- t(as.matrix(s))
  
  out <- FALSE
  if(length(s) == 1){ # if points are in 1D, check manually
    
    if(s[1] > min(A) & s[1] < max(A)) out <- TRUE
  }else{ # otherwise use inhulln from geometry package
    
    out <- inhulln(convhulln(A), s)
  }
  
  return(out)
}


get_ordered_hull <- function(convhull){
  
  n_vertices <- nrow(convhull)
  
  out <- c()
  current_vertex <- convhull[1, 1]
  which_edge <- 1
  while(length(out) < n_vertices){
    
    out <- append(out, current_vertex)
    current_edge <- convhull[which_edge, ]
    next_vertex <- current_edge[which(current_edge != current_vertex)]
    convhull <- convhull[!apply(convhull, 1, function(x) current_vertex %in% x), , drop = FALSE]
    current_vertex <- next_vertex
    which_edge <- which(apply(convhull, 1, function(x) current_vertex %in% x))
  }
  return(out)
}


wendell_prob <- function(m, p){
  
  ### calculate the corresponding Wendell's Theorem probability (Eq. 4) for a given m and p
  
  return(1 - sum(choose(m - 1, 0:(p - 2))) / 2^(m-1))
}

exponential_bound <- function(m, p){
  
  ### calculate exponential upper bound for hull probabilities (Eq. 5) for a given m and p
  
  return((1 - 2/(m + 1))^(p-1))
}
