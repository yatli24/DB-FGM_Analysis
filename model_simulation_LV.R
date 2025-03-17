# This script extends on the findings of the paper by running an increased
# sample size, curves, and changepoint to model the seasonality trend in Lake
# Victorial Basin. It stores the results into a folder of the current directory. 
# After running this script, the results are compiled into tables for 
# visualizations with the analyze_simulation_results.R script


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

# sample size, increased from paper's study
n = 100

# number of curves, increased from paper's study
p = 25

# set pii
pii_local = 1

# set number of bases needed for data generation
K_true = 5

# set basis as polynomial 
basis_type_true = 'polynomial'

# number of time points
T_data = 256

# set the true changepoint
changepoint_true = 160

# set continuous changepoint adjustment
continuous = TRUE    

# define standard noise error
sigma_epsilon_true = 0.05   

# simulate the data with function call
data_changepoint = simulate_data(p, K_true, n, T_data, pii_local, 
                                 sigma_epsilon_true, basis_type_true,
                                 continuous = TRUE, changepoint_true)

# save the simulated data as a Rdata object
folder_name = "Simulation_data"
file_name = paste("Simulation_data_new.Rdata", sep = "")
save(data_changepoint, file=paste(folder_name, '/', file_name, sep = ""))

# plot the new simulated data, save image
png("LV_simulation_time_series.png", width = 800, height = 800)
matplot(data_changepoint$Y[1,,], 
        type = 'l', 
        lty = 1, 
        col = 1:ncol(data_changepoint$Y[1,,]), 
        xlab = "Days After March", 
        ylab = "Simulated SST", 
        main = "Simulated Seasonlity Trend of Lake Victoria Basin SST")

dev.off()

# Run DBFGM
# Model parameters

# number of basis functions
K = 5
basis_type = 'spline'

# spike variance
v0 = 0.0004
h = 50^2;

# slab variance
v1 = v0 * h 
disp = TRUE

a_pi = 2; b_pi = 7
changepoint_interval = matrix(NA, 1, 2)
changepoint_interval[1,] = c(109, 149)  

# MCMC setup
nburn = 3000; nsave = 2000; 
nrep = 50
#save destination
folder_name = "Simulation_results_DBFGM"
dir.create(folder_name)

for (rep_ind in 1:nrep){ 
  print(rep_ind)
  
  # Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  
  set.seed(123 + rep_ind)
  # Run DBFGM
  mcmc_output = call_DBFGM(data_changepoint_rep,
                           K,
                           nburn, nsave,
                           v0, v1, a_pi, b_pi,
                           basis_type,
                           changepoint_interval,
                           disp)
  # Save results
  folder_name = "Simulation_results_DBFGM"
  file_name = paste("MCMC_output_DBFGM_rep", rep_ind, ".Rdata", sep = "")
  save(mcmc_output, file=paste(folder_name, '/', file_name, sep = ""))
  
  # Get model performances before and after the change point
  param_true = data_changepoint_rep$param_true
  for (s_i in 1:2){
    data_oneint = list()
    data_oneint$G_x_true = param_true$G_x_true[[s_i]]
    data_oneint$G_b_true = param_true$G_b_true[[s_i]]
    data_oneint$Omega_b_true = param_true$Omega_b_true[[s_i]]
    data_oneint$B_true = param_true$B_true[[s_i]]
    data_oneint$F_true = param_true$F_true
    data_oneint$K_true = param_true$K_true
    mcmc_output_oneint = mcmc_output[[s_i]]
    
    performance_graph = get_mcmc_perf_graph(mcmc_output_oneint, data_oneint, 
                                            K, p, standardize = FALSE, block_thresh = 0.5, disp = FALSE)
    cat("\nGraph estimation performance of state ", s_i, "\n")
    print_mcmc_results(performance_graph, data_oneint)
    
    file_name = paste("MCMC_performance_DBFGM_s", s_i, "_rep", rep_ind, ".Rdata", sep = "")
    save(performance_graph, file=paste(folder_name, '/', file_name, sep = ""))
  }
}

# after running the proposed model, the other models must be run via their respective scripts

# Run the partially separable functional graphical model (PSFGM)
file.edit('Helper_functions/Simulation_PSFGM_rep.R')

# Run BGGM
file.edit('Helper_functions/Simulations_BGGM_rep.R')

# Compare performance of models
file.edit('Helper_functions/analyze_simulation_results.R')
