# load necessary packages
library(BDgraph)
library(pheatmap)
library(fda)
library(coda)

# set a seed for reproducibility
set.seed(123)

# set the working directory to helper functions folder
setwd("Helper_functions")

# load necessary helper functions
source('Helper_functions_data_generation.R')
source('Helper_functions_MCMC.R')
source('performance_functions.R')
source('Call_functions.R')
source('MCMC_algorithms.R')

# reset directory
setwd("..")

# generate data, set up parameters

# sample size
n = 50

# number of curves
p = 15

# set pii
pii_local = 1

# set number of bases needed for data generation
K_true = 5

# set basis as polynomial 
basis_type_true = 'polynomial'

# number of time points
T_data = 256

# set the true changepoint
changepoint_true = 129

# set continuous changepoint adjustment
continuous = TRUE    

# define standard noise error
sigma_epsilon_true = 0.05   

# simulate the data with function call
data_changepoint = simulate_data(p, K_true, n, T_data, pii_local, 
                                 sigma_epsilon_true, basis_type_true,
                                 continuous = TRUE, changepoint_true)

# save the simulated data as a Rdata object
save(data_changepoint, file="Simulation_data/Simulation_data_original.Rdata")

# plot simulated data, Figure 2, save figures
# Simulated graph before changepoint
png("figures/graph_before_changepoint.png", width = 800, height = 800)
pheatmap(data_changepoint$param_true$G_x_true[[1]] + 0, 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"), 
         breaks = c(0, 0.5, 1))
dev.off()

# Simulated graph after changepoint
png("figures/graph_after_changepoint.png", width = 800, height = 800)
pheatmap(data_changepoint$param_true$G_x_true[[2]] + 0, 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"), 
         breaks = c(0, 0.5, 1))
dev.off()

# plot the simulated data time series
png("figures/simulated_time_series.png", width = 800, height = 600)
matplot(data_changepoint$Y[1,,], type = 'l', lty = 1, col = 1:ncol(data_changepoint$Y[1,,]))
dev.off()

