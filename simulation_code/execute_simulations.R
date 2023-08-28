### zrmiller@illinois.edu | zachary.miller@yale.edu ###

### Execute large-scale simulations

library(deSolve)
library(geometry)
library(tidyverse)
library(reshape2)

source("./simulation_functions.R")

###### Birth-death-immigration simulations (for Figs. 3 and S1) #####

# number of resources (to loop over) 
p_vec <- c(2, 4, 8, 16)

m_init <- 2 # initial number of consumer species 
sigma <- 0 # standard deviation of tradoff noise (0 for perfect trade-off; 0.0005 for Fig S1)
mutation_rate <- 0.001 # probability that a new species enter the system at each step
S <- 5000 # total resource supply = total number of consumer individuals
max_steps <- 10^7 # number of total time steps
n_reps <- 50 # number of independent realizations of the dynamics (each with a fixed s vector sampled at random)

# initialize tibble for results
results <- tibble(p = numeric(), rep = numeric(), timepoint = numeric(), abuns = numeric())

# for paired neutral simulations
results_neutral <- tibble(p = numeric(), rep = numeric(), timepoint = numeric(), abuns = numeric())

for(i in 1:n_reps){
  for(p in p_vec){
    
    print(paste("p:", p, "replicate:", i)) # print progress
    
    # generate parameters
    A_exact <- generate_tradeoff_matrix(p, m_init) # sample matrix of resource utilization rates
    A <- A_exact + matrix(rnorm(p * m_init, 0, sigma), nrow = m_init, ncol = p) # add noise to entries
    
    s <- sample_from_simplex(p) # sample resource supply uniformly from simplex (neglect scaling by S, as this does not change dynamics)
    
    n_init <- S * rep(1/m_init, m_init) # initialize all consumers with equal abundances 
    
    pars <- list(A = A, s = s, sigma = sigma) # combine all parameters
    
    # simulate BDI model
    out_stochastic <- stochastic_birth_death_immigration_dynamics(state = n_init, pars = pars,
                                                                  mutation_rate = mutation_rate, max_steps = max_steps)
    
    results <- results %>% 
      add_row(out_stochastic %>% 
                mutate(p = p, rep = i)
      ) # add results for this replicate to full set of results
  }
  
  # simulate one set of replicates with a purely neutral model for comparison
  out_neutral <- stochastic_birth_death_immigration_dynamics(state = n_init, pars = "neutral", 
                                                             mutation_rate = mutation_rate, max_steps = max_steps)
  
  results_neutral <- results_neutral %>% 
    add_row(out_neutral %>% 
              mutate(p = 1, rep = i)
    ) # add results for this replicate to full set of results
}

results_all <- rbind(results_neutral, results) # combine neutral and non-neutral model results

# write output to csv
write_csv(x = results_all, 
          file = paste0("../simulation_results/BDI_nu_", mutation_rate, "_S_", S, "_reps_", n_reps, "_sigma_", sigma, ".csv"))


###### Birth-death simulations (for Figs. 4 and 5) #####

# number of resources (to loop over)
p_vec <- 2:9 

m <- 10 # number of consumer species
sigma <- 0 # standard deviation of tradoff noise
S <- 5000 # total resource supply

n_reps <- 1000 # how many random parameterizations for p

# initialize matrix for simulation results (2 rows for each p x rep combination, due to 2 resource supply scenarios)
results <- matrix(nrow = 2 * length(p_vec) * n_reps,
                  ncol = 7)
counter <- 1 # counter for row in results

# for paired neutral simulations
results_neutral <- matrix(nrow = 2 * n_reps,
                          ncol = 7)
neutral_counter <- 1 # counter for row in results_neutral

for(i in 1:n_reps){
  for(p in p_vec){
    for(s_type in 1:2){
      
      print(paste("p:", p, "replicate:", i)) # print progress
      
      # generate parameters
      A_exact <- generate_tradeoff_matrix(p, m) # sample matrix of resource utilization rates
      A <- A_exact + matrix(rnorm(p * m_init, 0, sigma), nrow = m_init, ncol = p) # add noise to entries
      
      if(s_type == 1){
        s <- rep(1/p, p) # equal supply of all resources (neglect scaling by S, as this does not change dynamics)
      }else{
        s <- sample_from_simplex(p) # sample resource supply uniformly from simplex (neglect scaling by S, as this does not change dynamics)
      }
      
      n_init <- S * rep(1/m, m) # initialize all consumers with equal abundances 
      
      pars <- list(A = A, s = s)
      
      inhull <- test_inhull(A_exact, s, simplex = TRUE) # check if the resource supply vector is in the convex hull of consumer metabolic strategies
      
      # simulate BD model
      out_stochastic <- stochastic_birth_death_dynamics(state = n_init, pars = pars,
                                                        until = 1) # simulate until first extinction occurs
      
      results[counter, ] <- c(p, m, sigma, i, inhull, s_type,
                              out_stochastic[1] / S) # divide time by S to normalize by total system size
      counter <- counter + 1
    }
  }
  
  # simulate one set of replicates with a purely neutral model for comparison
  
  out_neutral <- stochastic_birth_death_dynamics(state = n_init, pars = "neutral", 
                                                 until = 1) # simulate until first extinction occurs
  
  results_neutral[neutral_counter, ] <- c(1, m, NA, i, NA, 1, out_neutral[1] / S) # s_type = 1
  results_neutral[neutral_counter + 1, ] <- c(1, m, NA, i, NA, 2, out_neutral[1] / S) # s_type = 2
  # NOTE: s_type has no effect on neutral model, but both types are needed for plotting -- to account for this, save two copies of neutral results
  neutral_counter <- neutral_counter + 2
}

# organize simulation results for plotting
results <- results %>% as_tibble() %>% 
  setNames(c("p", "m", "sigma", "rep", "inhull", "s_type", "first_extinction"))
results_neutral <- results_neutral %>% as_tibble() %>%
  setNames(c("p", "m", "sigma", "rep", "inhull", "s_type", "first_extinction"))

# write output to csv
write_csv(x = results, 
          file = paste0("../simulation_results/BD_S_", S, "_reps_", n_reps, "_sigma_", sigma, ".csv"))
write_csv(x = results, 
          file = paste0("../simulation_results/BD_S_", S, "_reps_", n_reps, "_neutral.csv"))


##### Compute convex hull probabilities (for Fig. 6) #####

# number of consumers (to loop over)
m_vec <- seq(10, 25, by = 5)

n_samples <- 1000
s_type <- "equal" # options are "equal" or "random"

# initialize matrix for results
hull_probabilities <- matrix(nrow = length(m_vec), 
                             ncol = max(m_vec) - 2)
m_counter <- 0 # to keep track of m value
for(m in m_vec){
  
  m_counter <- m_counter + 1
  
  p_counter <- 0 # to keep track of p value
  for(p in 2:(m-1)){
    
    p_counter <- p_counter + 1
    print(paste("m:", m, "p:", p)) # print progress
    
    hull_counter <- 0 # to count number of parameters in hull
    for(i in 1:n_samples){
      
      A <- generate_tradeoff_matrix(p, m) # generate a random resource utilization matrix
      
      # now sample a random resource supply vector
      if(identical(s_type, "equal")){
        s <- rep(1/p, p) # equal supply of all resources
      }else{
        if(identical(s_type, "random")){
          s <- sample_from_simplex(p) # sample resource supply uniformly from simplex 
        }else{
          print("Invalid s_type")
          break
        }
      }
      
      inhull <- test_inhull(A, s, simplex = TRUE) # check if s is in the hull of A
      if(inhull) hull_counter <- hull_counter + 1
    }
    hull_probabilities[m_counter, p_counter] <- hull_counter
  }
}

hull_probabilities <- hull_probabilities %>% melt() # convert results to tidy form
colnames(hull_probabilities) <- c("m", "p", "probability")
hull_probabilities <- hull_probabilities %>% 
  mutate(m = m_vec[m], # map m and p counters to correct values
         p = p + 1, 
         probability = probability / n_samples) %>% # divide in-hull count by total samples to get a probability
  filter(!is.na(probability)) # remove unused m / p combinations

# select appropriate (vectorized) theoretical prediction for comparison with simulation results
if(identical(s_type, "equal")){
  theory_prediction_func <- Vectorize(wendell_prob)
}else{
  if(identical(s_type, "random")){
    theory_prediction_func <- Vectorize(exponential_bound)
  }else{
    print("Invalid s_type")
    break
  }
}

hull_probabilities <- hull_probabilities %>% mutate(theory_prediction = theory_prediction_func(m, p)) # add theoretical prediction

write_csv(x = hull_probabilities, 
          file = paste0("../simulation_results/hull_probabilities_", s_type, "_samples_", n_samples, ".csv"))
