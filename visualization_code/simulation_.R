### zrmiller@illinois.edu | zachary.miller@yale.edu ###

### Generate plots shown in the paper

library(deSolve)
library(Ternary)
library(geometry)
library(tidyverse)
library(reshape2)
library(scales)

source("./simulation_functions.R")

##### Fig 1 : convex hull condition and corresponding dynamics #####

### First, simulate dynamics with and without convex hull condition

seed_num <- 15
set.seed(seed_num)

p <- 3 # number of resources
m <- 4 # number of species

# total resource supply
S <- 2000

# generate resource utilization matrix
A <- generate_tradeoff_matrix(p, m)

# For panel b, use a resource supply outside the convex hull (using seed 15)
s <- S * sample_from_simplex(p) # sample resource supply uniformly from simplex 
sb <- s # a copy of this resource supply for later plotting

# initial conditions
n_init <-  S * c(0.3, 0.4, 0.2, 0.1) 

# times to simulate
times <- seq(0, 50, by = 1/S)

reduced_pars <- list(A = A, s = s)

# simulate simplified ODE model
out_reduced <- ode(y = n_init,
                   times = times,
                   func = consumer_dynamics,
                   parms = reduced_pars,
                   method = "ode45")

# simulate stochastic BD model with same parameters
out_stochastic <- stochastic_birth_death_dynamics(time = times, state = n_init, pars = reduced_pars)

# organize results
consumers_reduced <- out_reduced %>% as.data.frame() %>% 
  mutate(method = "reduced", type = "consumer") %>%
  pivot_longer(cols = -c("time", "method", "type"))
consumers_stochastic <- out_stochastic %>% as.data.frame() %>% 
  mutate(method = "stochastic", type = "consumer") %>%
  pivot_longer(cols = -c("time", "method", "type"))

all_tidy <- rbind(consumers_reduced, consumers_stochastic)

# plot consumer abundances for both models (ODE + stochastic BD)
p_1b <- all_tidy %>% 
  filter(value > 0) %>%
  ggplot() + 
  aes(x = time, y = value, 
      group = interaction(name, method), # name is species, method is model type
      color = name, alpha = method) + 
  geom_line() + 
  scale_color_viridis_d(guide = "none") + 
  theme_bw() + 
  scale_alpha_manual(values = c(1, 0.5), guide = "none") + 
  scale_y_log10(limits = c(1, S)) + 
  xlab("Time") + 
  ylab("Abundance")

# next, for panel c, use a resource supply inside the convex hull (using seed 15)

s <- S * c(0.25, 0.45, 0.3)
sc <- s # a copy of this resource supply for later plotting

# initial conditions (use same as for panel b)
#n_init <- S * c(0.3, 0.4, 0.2, 0.1)

# times to simulate (use same as for panel b)
#times <- seq(0, 50, by = 1/S)

reduced_pars <- list(A = A, s = s)

# simulate simplified ODE model
out_reduced <- ode(y = n_init,
                   times = times,
                   func = consumer_dynamics,
                   parms = reduced_pars,
                   method = "ode45")

# simulate stochastic BD model with same parameters
out_stochastic <- stochastic_birth_death_dynamics(time = times, state = n_init, pars = reduced_pars)

# organize results
consumers_reduced <- out_reduced %>% as.data.frame() %>% 
  mutate(method = "reduced", type = "consumer") %>%
  pivot_longer(cols = -c("time", "method", "type"))
consumers_stochastic <- out_stochastic %>% as.data.frame() %>% 
  mutate(method = "stochastic", type = "consumer") %>%
  pivot_longer(cols = -c("time", "method", "type"))

all_tidy <- rbind(consumers_reduced, consumers_stochastic)

# plot consumer abundances for both models (ODE + stochastic BD)
p_1c <- all_tidy %>% 
  filter(value > 0) %>%
  ggplot() + 
  aes(x = time, y = value, 
      group = interaction(name, method), 
      color = name, alpha = method) + 
  geom_line() + 
  scale_color_viridis_d(guide = "none") + 
  theme_bw() + 
  scale_alpha_manual(values = c(1, 0.5), guide = "none") + 
  scale_y_log10(limits = c(1, S)) + 
  xlab("Time") + 
  ylab("Abundance")

### For panel a, plot metabolic simplex with consumer strategies and both resource supply points

TernaryPlot(grid.lines = 0) # blank ternary plot
ordered_hull <- A[get_ordered_hull(convhulln(A[, -1])), ] # get consumer strategies in a format for plotting
TernaryPolygon(ordered_hull, col = "grey", border = "grey") # shade in convex hull
AddToTernary(points, A, col = viridis(m), pch = 20) # add consumer strategies as colored points 
AddToTernary(text, sb, "b")
AddToTernary(text, sc, "c")


##### Fig 2 : emergent neutrality #####

seed_num <- 10
set.seed(seed_num)

p <- 2 # number of resources
m <- 3 # number of species

# generate resource utilization matrix
A <- generate_tradeoff_matrix(p, m)

# total resource supply
S <- 2000
s <- S * rep(1/p, p) # equal supply of all resources

# initial conditions
n_init <- S * c(0.7, 0.1, 0.2)

# times to simulate
times <- seq(0, 200, by = 1/S)

reduced_pars <- list(A = A, s = s)

# simulate simplified ODE model
out_reduced <- ode(y = n_init,
                   times = times,
                   func = consumer_dynamics,
                   parms = reduced_pars,
                   method = "ode45")

# simulate stochastic BD model with same parameters
out_stochastic <- stochastic_birth_death_dynamics(time = times, state = n_init, pars = reduced_pars)

# organize results
consumers_reduced <- out_reduced %>% as.data.frame() %>% 
  mutate(method = "reduced", type = "consumer") %>%
  pivot_longer(cols = -c("time", "method", "type"))
consumers_stochastic <- out_stochastic %>% as.data.frame() %>% 
  mutate(method = "stochastic", type = "consumer") %>%
  pivot_longer(cols = -c("time", "method", "type"))

all_tidy <- rbind(consumers_reduced, consumers_stochastic)

# plot consumer abundances for both models (ODE + stochastic BD)
p_2a <- all_tidy %>% filter(time > 0.5) %>%
  ggplot() + 
  aes(x = time, y = value, group = interaction(name, method), color = name, alpha = method) + 
  geom_line() + 
  scale_color_viridis_d() + 
  scale_alpha_manual(values = c(1, 0.5)) + 
  geom_vline(xintercept = 4, linetype = "dotted") + 
  theme_bw() + 
  xlab("Time") + 
  ylab("Abundance") + 
  scale_x_log10() + 
  theme(legend.position = "none")

### plot same dynamics in the simplex of consumer abundances

# generate the neutral manifold for these parameters
length.out <- 100
neutral_manifold <- matrix(nrow = length.out, ncol = 3)
for(i in 1:length.out){
  n1 <- S * (i / length.out)
  n2 <- ((1 - A[3, 1]) * s[1] - A[3, 1] * s[2] - (A[1, 1] - A[3, 1]) * n1) / (A[2, 1] - A[3, 1])
  n3 <- S - n1 - n2
  
  neutral_manifold[i, ] <- c(n1, n2, n3)
}
feasible_points <- apply(neutral_manifold, 1, function(x) all(x > 0)) # keep only points that are feasible
neutral_manifold <- neutral_manifold[feasible_points, ]

# plot in the simplex
TernaryPlot(grid.lines = 0, padding = 0)
AddToTernary(lines, out_reduced[, -1]) # plot ODE dynamics (reduced dimension)
AddToTernary(lines, neutral_manifold, lty = "dashed") # show neutral manifold
AddToTernary(lines, # plot stochastic dynamics (reduced dimension)
             out_stochastic[seq(from = 1, to = nrow(out_stochastic), length.out = 1000), -1],
             col = rgb(red = 0, green = 0 , blue = 1, alpha = 0.5))
AddToTernary(points, out_stochastic[4 * S, -1]) # annotate approximate transition to emergent neutrality
AddToTernary(points, n_init, pch = 4) # annotate initial conditions


##### Fig 3 : BDI simulations #####

BDI_dat <- read_csv("../simulation_results/BDI_nu_0.001_S_5000_sigma_0.csv")
BDI_dat <- BDI_dat %>% filter(timepoint > 50) # use only late times to ensure convergence to steady state

# theoretical formula for expected number of species in neutral model
ewens_formula <- function(n, S, nu){
  theta <- (S - 1) * (nu / (1 - nu))
  return(exp(log(theta) + lgamma(S + 1) + lgamma(S + theta - n) - (log(n) + lgamma(S + 1 - n) + lgamma(S + theta))))
}

# compute rank abundance curves for plotting
prepare_rank_curves <- function(x){
  
  out <- x %>%
    # for each p, timepoint, and replicate, rank species by abundance
    group_by(timepoint, rep, p) %>% mutate(rank = rank(-abuns, ties.method = "first")) %>% ungroup() %>% 
    mutate(n_times = max(timepoint) - min(timepoint),
           n_reps = max(rep)) %>%
    group_by(p, rank) %>% 
    summarize(mean_abun = sum(abuns) / (n_reps * n_times))
  
  return(out)
}

p_3a <- BDI_dat %>%
  filter(p > 1) %>% # exclude neutral model for now and plot separately
  prepare_rank_curves %>%
  ggplot() + 
  aes(x = rank, y = mean_abun) + 
  geom_line(aes(color = as.factor(p)), show.legend = FALSE) +
  geom_tile(data = data.frame(rank = rep(Inf, 4), mean_abun = rep(Inf, 4), p = c(2, 4, 8, 16)), # plot dummy points to produce color bar
            aes(fill = as.factor(p))) +
  scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x))) +
  xlim(0, 60) + 
  scale_fill_viridis_d(name = "Number of \n resources") + 
  scale_color_viridis_d(name = "Number of \n resources") + 
  theme_bw() + 
  theme(legend.position = c(0.8, 0.6),
        legend.title = element_text(size = 10),
        legend.margin = margin(c(0,0.1,0.1,0.1)),
        legend.key.size = unit(0.4, "cm")) + 
  xlab("Rank") + ylab("Mean abundance") 

p_3b <- BDI_dat %>%
  filter(p > 1) %>% # exclude neutral model for now and plot separately
  filter(abuns > 1) %>% # exclude singletons
  group_by(p, timepoint, rep) %>% count() %>% # get richness for each p, timepoint, and replicate
  ggplot() + 
  aes(x = as.factor(p), y = n, color = as.factor(p)) + 
  geom_hline(yintercept = -S * mutation_rate * log(mutation_rate) - ewens_formula(1, S, mutation_rate), # plot theoretical richness (adjusted to remove singletons)
             linetype = "dashed") + 
  geom_boxplot(outlier.alpha = 0.25) +
  scale_color_viridis_d() + 
  xlab("Number of resources") + 
  ylab("Number of species") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw() + 
  theme(legend.position = "none")

ggarrange(plotlist = list(p_3a, p_3b), labels = "auto", widths = c(1.5, 1))


##### Fig 4 : BD simulations; marginal CDFs #####

BD_dat <- read_csv("../simulation_results/BD_S_5000_reps_3000_sigma_0.csv")

# convenience function to compute CDF
my_ecdf <- function(x) ifelse(x == 0, 0, ecdf(x)(x))

p_4 <- BD_dat %>% 
  group_by(s_type, p) %>% mutate(cdf = my_ecdf(first_extinction)) %>% ungroup(p) %>% # compute CDFs
  ggplot() + aes(x = first_extinction, y = 1 - cdf, group = as.factor(p), fill = as.factor(p)) + 
  geom_line(aes(color = as.factor(p)), show.legend = FALSE) + 
  geom_tile(data = data.frame(first_extinction = rep(Inf, 8), cdf = rep(-Inf, 8), p = 2:9)) + # plot dummy points to produce color bar
  facet_grid(.~s_type, scales = "free",
             labeller = as_labeller(c(`1` = "Equal", `2` = "Random"))) + # label resource supply scenarios
  scale_fill_viridis_d(name = NULL) + 
  scale_color_viridis_d() + 
  theme_bw() + 
  theme(legend.position = c(0.96, 0.69), 
        legend.margin = margin(c(0,0.1,0.1,0.1)),
        legend.key.size = unit(0.4, "cm")) + 
  scale_x_log10() +
  ylab("Fraction of communities") + xlab("Time to extinction") 


##### Fig 5 : BD simulations; conditional and marginal times #####

BD_dat <- read_csv("../simulation_results/BD_S_5000_reps_3000_sigma_0.csv")
BD_neutral_dat <- read_csv("../simulation_results/BD_S_5000_reps_1000_neutral.csv")

# info for manual legend
legend_df <- data.frame(p = rep(6.2, 3), first_extinction = c(25, 30, 36), 
                        s_type = rep(1, 3), inhull = c("Conditional: In", "Conditional: Out", "Marginal"))

p_5 <- BD_dat %>% 
  mutate(inhull = ifelse(inhull, "Conditional: In", "Conditional: Out")) %>% # map inhull to a more informative label
  add_row(BD_dat %>% mutate(inhull = "Marginal")) %>% # create an unconditioned (marginal) copy of results
  group_by(p, s_type, inhull) %>% mutate(n = n()) %>% ungroup() %>% # count observations in each category (in hull vs. out)
  ggplot() + aes(x = p, first_extinction) + 
  geom_hline(yintercept = quantile(BD_neutral_dat$first_extinction, 0.5), linetype = "dashed") + # add median time to first extinction in neutral model
  annotate(geom = "rect", # add interquartile range for neutral model
           ymin = quantile(BD_neutral_dat$first_extinction, 0.25),
           ymax = quantile(BD_neutral_dat$first_extinction, 0.75),
           xmin = -Inf,
           xmax = +Inf,
           fill = "lightgray", alpha = 0.5, size = 0) +
  stat_summary(geom = "errorbar", aes(group = inhull, color = inhull), # plot interquartile range for conditional distributions
               fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75), 
               width = 0,
               position = position_dodge(width = 0.25),
               show.legend = FALSE) + 
  stat_summary(geom = "point", # plot medians for conditional distributions
               aes(group = inhull, color = inhull, size = n),
               fun = function(x) quantile(x, 0.5), 
               position = position_dodge(width = 0.25), shape = 19) + 
  stat_summary(geom = "errorbar", aes(alpha = inhull), # plot interquartile range for marginal distributions
               fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75), 
               width = 0,
               position = position_dodge(width = 0.25),
               color = "#440154FF",
               show.legend = FALSE) + 
  geom_point(data = legend_df, aes(color = inhull)) + # create manual legend
  geom_label(data = legend_df, aes(x = p + 0.25, label = inhull), # create manual legend
             hjust = "left", size = 3, label.size = 0) +
  facet_wrap(.~s_type, scales = "free",
             labeller = as_labeller(c(`1` = "Equal", `2` = "Random"))) + # label resource supply scenarios
  scale_color_viridis_d(direction = -1, guide = "none") + 
  scale_alpha_manual(values = c(0, 0, 1), guide = "none") + 
  scale_size(range = c(0.5, 3), guide = "none") + 
  theme_bw() +
  scale_y_log10() + 
  xlab("Number of resources") + ylab("Time to extinction")

show(p)