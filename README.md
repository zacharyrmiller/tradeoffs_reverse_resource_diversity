Code to re-create figures and analysis from "Metabolic trade-offs can reverse the resource-diversity relationship" by Zachary R. Miller and James P. O'Dwyer (The American Naturalist 2024).

Email ZR Miller with any questions (zachary.miller@yale.edu).

Short abstract:
This study analyzes a mathematical model for consumer-resource dynamics with metabolic trade-offs and shows that trade-offs can reverse the expected relationship between resource diversity and consumer diversity. This reversal is a result of a generic transition between neutral-like and niche-like dynamics, driven by resource diversity.

Written and maintained by ZR Miller.

Study results and figures were produced with R version 4.3.2.

Required packages:

ggpubr_0.6.0      viridis_0.6.5     viridisLite_0.4.2 scales_1.3.0     
Ternary_2.3.1     reshape2_1.4.4    lubridate_1.9.3   forcats_1.0.0    
stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.5      
tidyr_1.3.1       tibble_3.2.1      ggplot2_3.4.4     tidyverse_2.0.0  
geometry_0.4.7    deSolve_1.40 

This repository includes three folders:

- simulation_code -- contains functions (simulation_functions.R) and scripts to execute (execute_simulations.R) all simulations reported in the paper
  
- visualization code -- contains code used to generate all figures in the paper (figures.R)
  
- simulation_results -- all simulation results used in the paper, in csv format
    - Files of the format BDI_nu_X_S_XX_reps_XXX_XXXX.csv record results for the birth-death-immigration model with nu (mutation/immigration rate) equal to X, total resource inflow (equivalently, system size) equal to XX, and sigma (noise around the trade-off) equal to XXXX (or no tradeoff, as indicated), for XXX replicate simulations. Data are recorded in long (tidy) form, with the variables (columns): p (number of resources), rep (replicate label), timepoint (timestep in the simulation), and abuns (abundance of a consumer species in the simulation). For each combination of p, rep, and timepoint, distinct values of abuns indicate distinct species (species identities not tracked through time).
    - Files of the format BD_S_X_reps_XX_XXX.csv record results for the birth-death model with total resource inflow (equivalently, system size) equal to X and sigma (noise around the trade-off) equal to XXX (or no tradeoff, as indicated), for XX replicate simulations. Variables (columns) are: p (number of resources), m (initial number of consumer species), sigma (noise around the tradeoff, if relevant), rep (replicate label), inhull (logical variable indicating whether (1) or not (0) the resource supply vector is within the convex hull of consumer metabolic strategies), s_type (distribution of resource supply, with 1 = equal supply and 2 = random supply), and first extinction (time until at least one consumer species went extinct -- arbitrary time units)
