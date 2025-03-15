# After running all 3 models on the simulated data,
# this script compiles the results and saves them into RData tables

library(BDgraph)
library(pheatmap)
library(fda)
library(coda)
library(knitr)

rm(list = ls()); 
setwd("C:/Users/ivani/OneDrive/Documents/Offline R project/DBFGM-main/Helper_functions")
source('performance_functions.R')
setwd("..")

# Load results for proposed model
nrep = 50
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  for (s_i in 1:2){
    folder_name = "Simulation_results_DBFGM"
    file_name = paste("MCMC_performance_DBFGM_s", s_i,"_rep", rep_ind, ".Rdata", sep = "")
    load(file=paste(folder_name, '/', file_name, sep = ""))
    tpr[rep_ind, s_i] = performance_graph$block_perf$tpr_block
    fpr[rep_ind, s_i] = performance_graph$block_perf$fpr_block
    mcc[rep_ind, s_i] = performance_graph$block_perf$mcc_block
    
  }
}

# Compute statistics
tpr_mean <- apply(tpr, 2, mean)
tpr_sd <- apply(tpr, 2, sd)
fpr_mean <- apply(fpr, 2, mean)
fpr_sd <- apply(fpr, 2, sd)
mcc_mean <- apply(mcc, 2, mean)
mcc_sd <- apply(mcc, 2, sd)

# Create a table
dbfgm_table <- data.frame(
  Metric = rep(c("TPR", "FPR", "MCC"), each = 4),
  Statistic = rep(c("Mean", "SD"), times = 6),
  Scenario = rep(c("s1", "s2"), each = 6),
  Value = c(tpr_mean, tpr_sd, fpr_mean, fpr_sd, mcc_mean, mcc_sd)
)

# save the results dataframe
save(dbfgm_table, file='table_dbfgm_results.RData')

# Create the kable object
kable_dbfgm <- kable(dbfgm_table, caption = "Performance Metrics Summary of DBFGM")

# Save the kable object as an RData file
save(kable_dbfgm, file = "kable_dbfgm_results.RData")


# Load results PSFGM
folder_name = "Simulation_results_PSFGM/"
nrep = 50
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  for (s_i in 1:2){
    file_name = paste("Output_PSFGM_rep", rep_ind, ".Rdata", sep = "")
    load(file=paste(folder_name, '/', file_name, sep = ""))
    tpr[rep_ind, s_i] = psfgm_output$performance[[s_i]][1]
    fpr[rep_ind, s_i] =  psfgm_output$performance[[s_i]][2]
    mcc[rep_ind, s_i] =  psfgm_output$performance[[s_i]][3]
    
  }
}

# Compute statistics
tpr_mean <- apply(tpr, 2, mean)
tpr_sd <- apply(tpr, 2, sd)
fpr_mean <- apply(fpr, 2, mean)
fpr_sd <- apply(fpr, 2, sd)
mcc_mean <- apply(mcc, 2, mean)
mcc_sd <- apply(mcc, 2, sd)

# Create a table
psfgm_table <- data.frame(
  Metric = rep(c("TPR", "FPR", "MCC"), each = 4),
  Statistic = rep(c("Mean", "SD"), times = 6),
  Scenario = rep(c("s1", "s2"), each = 6),
  Value = c(tpr_mean, tpr_sd, fpr_mean, fpr_sd, mcc_mean, mcc_sd)
)

# save the results dataframe
save(psfgm_table, file='table_psfgm_results.RData')

# Create the kable object
kable_psfgm <- kable(psfgm_table, caption = "Performance Metrics Summary of PSFGM")

# Save the kable object as an RData file
save(kable_psfgm, file = "kable_psfgm_results.RData")

# BGGM -----------------
nrep = 50
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  folder_name = "Simulation_results_BGGM"
  file_name = paste("MCMC_performance_BGGM_rep", rep_ind, ".Rdata", sep = "")
  load(file=paste(folder_name, '/', file_name, sep = ""))

  for (s_i in 1:2){
    tpr[rep_ind, s_i] = performance_graph[[s_i]][1]
    fpr[rep_ind, s_i] =  performance_graph[[s_i]][2]
    mcc[rep_ind, s_i] =  performance_graph[[s_i]][3]
  }
  
}

# Compute statistics
tpr_mean <- apply(tpr, 2, mean)
tpr_sd <- apply(tpr, 2, sd)
fpr_mean <- apply(fpr, 2, mean)
fpr_sd <- apply(fpr, 2, sd)
mcc_mean <- apply(mcc, 2, mean)
mcc_sd <- apply(mcc, 2, sd)

# Create a table
bggm_table <- data.frame(
  Metric = rep(c("TPR", "FPR", "MCC"), each = 4),
  Statistic = rep(c("Mean", "SD"), times = 6),
  Scenario = rep(c("s1", "s2"), each = 6),
  Value = c(tpr_mean, tpr_sd, fpr_mean, fpr_sd, mcc_mean, mcc_sd)
)

# save the results dataframe
save(bggm_table, file='table_bggm_results.RData')

# Create the kable object
kable_bggm <- kable(bggm_table, caption = "Performance Metrics Summary of BGGM")

# Save the kable object as an RData file
save(kable_bggm, file = "kable_bggm_results.RData")



