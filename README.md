Code to re-create figures and analysis from "Metabolic trade-offs can reverse the resource-diversity relationship" by Zachary R. Miller and James P. O'Dwyer.

Email ZR Miller with any questions (zrmiller@illinois.edu OR zachary.miller@yale.edu).

This study analyzes a mathematical model for consumer-resource dynamics with metabolic trade-offs and shows that trade-offs can reverse the expected relationship between resource diversity and consumer diversity. This reversal is a result of a generic transition between neutral-like and niche-like dynamics, driven by resource diversity.

Written and maintained by ZR Miller.

Required packages:

deSolve, Ternary, geometry, tidyverse, reshape2, scales

This repository includes three folders:

- simulation_code -- contains functions (simulation_functions.R) and scripts to execute (execute_simulations.R) all simulations reported in the paper
- simulation_results -- all simulation results used in the paper, in csv format
- visualization code -- contains code used to generate all figures in the paper (figures.R)
