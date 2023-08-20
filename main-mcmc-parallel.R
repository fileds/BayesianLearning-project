# Loading packages
library(tidyverse)
library(mybnlearn)
library(foreach)
library(doParallel)
library(beepr)

rm(list = ls())

# Sourcing functions
source("R/dag_mcmc.R")
source("R/utils.R")

set.seed(42)
verbose = TRUE
sound = TRUE

# Loading the ASIA data from bnlearn.
data(asia)
# Asia true network
true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
# Sampling a subset of the data.
sample_size <- 100
#sample_size <- 1000
# Sampling data
df <- asia[sample(1:nrow(asia), sample_size), ]

# NUmber of correct edges
n_correct <- 6
 
# Prior imaginary sample size: 
# log(posterior) \propto loglik(D | G) + \prior_iss * P(G | \Xi)
# where D is the data, G is the DAG (structure), \prior_iss the imaginary 
# sample size for the prior, \Xi is the prior knowledge.
prior_isss <- c(1, 5, 10)
#prior_isss <- c(1, 10, 100)

# Causal Discovery parameters
# Implementation of score function.
score_fn <- mybnlearn::score

# Prior specification
# Complete prior specification with lower probability of inclusion for edges not
# in true and higher for edges in true.
prior <- prior_from_dag(true, prob = c(0.6, 0.2))

# MCMC
# Parameters
N <- 5000 # Number of MCMC steps

burnin <- ceiling(0.2 * N) # Burn in for the MCMC simulation.

# Simulation study
n_chains <- 8

# Setting up cluster.
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

t_start <- Sys.time()
result <- foreach(chain = 1:n_chains, .combine = rbind, .packages = c("tibble", "mybnlearn", "digest")) %dopar%
{
  cat("Processing chain ", chain, "\n")
  
  # Starting from a DAG with half of the edges correct
  start <- randomize_dag_edges(true, n_correct)
  #start <- mmhc(df, maximize.args = list(score = "bde"))
  start_score <- score(start, df, type = "bde")
  
  local_result <- tibble(
    step = 0,
    hash = NA,
    dag = list(start),
    score = start_score,
    prob = NA,
    accepted = NA,
    phase = "Start",
    prior = NA,
    iss = NA,
    chain = chain
  )
  for (i in 1:length(prior_isss))
  {
    iss <- prior_isss[i]
    if (verbose) cat(sprintf("iss = %d", iss))
    
    # MCMC with informative prior
    posterior <- dag_mcmc(
      df, 
      N, 
      score_fn,
      start = start,
      type = "bde", 
      prior = "cs", 
      beta = prior,
      prior_iss = iss
    )
    posterior$phase <- c(rep("Burn", burnin), rep("Sampling", N - burnin)) 
    posterior$prior <- rep("CS", nrow(posterior)) 
    posterior$iss <- rep(iss, nrow(posterior)) 
    posterior$chain <- rep(chain, nrow(posterior)) 
    
    local_result <- rbind(local_result, posterior)
  }
  # MCMC with uniform prior
  posterior <- dag_mcmc(
    df, 
    N, 
    score_fn,
    start = start,
    type = "bde", 
    prior = "uniform", 
  )
  posterior$phase <- c(rep("Burn", burnin), rep("Sampling", N - burnin)) 
  posterior$prior <- rep("Uniform", nrow(posterior)) 
  posterior$iss <- rep(0, nrow(posterior)) 
  posterior$chain <- rep(chain, nrow(posterior)) 
  
  local_result <- rbind(local_result, posterior)
}
t_end <- Sys.time()
cat(paste0("Time elapsed: ", format_elapsed_time(t_start, t_end), "\n"))

# Closing cluster
stopCluster(cl)

# Generate a timestamp for the filename
timestamp <- format(Sys.time(), format = "%Y-%m-%d_%H-%M-%S")
# Create the filename
filename <- paste("simulations/mcmc/mcmc_", sample_size, "_", timestamp, ".rds", sep = "")
# Save the object as an .rds file with the timestamped filename
saveRDS(result, file = filename)
if(sound) beep(sound = "fanfare")