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
#sample_size <- 200
#sample_size <- 500
#sample_size <- 1000
#sample_size <- 5000

# Scoring DAGs with CS prior
# Defining prior as true network
# Complete prior specification with lower probability of inclusion for edges not
# in true and higher for edges in true.
prior <- prior_from_dag(true, prob = c(0.6, 0.2))
 
# Prior imaginary sample size: 
# log(posterior) \propto loglik(D | G) + \beta_iss * P(G | \Xi)
# where D is the data, G is the DAG (structure), \beta_iss the imaginary 
# sample size for the prior, \Xi is the prior knowledge.
prior_isss <- c(1, 3, 5, 50)
#prior_isss <- c(1, 5, 10, 100)
#prior_isss <- c(1, 25, 250)
#prior_isss <- c(1, 50, 500)
#prior_isss <- c(1, 250, 500, 1000, 2500)

# MCMC
# Parameters
N <- 5000 # Number of MCMC steps
burnin <- ceiling(0.1 * N) # Burn in for the MCMC simulation.

# Simulation study
n_experiments <- 10

# Setting up cluster.
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

t_start <- Sys.time()
result <- foreach(experiment = 1:n_experiments, .combine = rbind, .packages = c("tidyverse", "mybnlearn", "digest")) %dopar%
{
  df <- asia[sample(1:nrow(asia), sample_size), ]
  
  # Causal Discovery parameters
  # Identifying starting DAG using Hill-Climbing algorithm with uniform BDe 
  # score.
  #start <- randomize_dag_edges(true, floor(nrow(arcs(true)) / 2))
  #start <- hc(df, score = "bde")
  start <- mmhc(df, maximize.args = list(score = "bde"))
  
  # MCMC with uniform prior
  posterior <- dag_mcmc(
    df, 
    N, 
    mybnlearn::score,
    start = start,
    type = "bde", 
    prior = "uniform", 
  )
  posterior <-
    posterior %>% 
    mutate(phase = c(rep("Burn", burnin), rep("Sampling", N - burnin)))
  
  # Count individual DAGs
  posterior <- posterior %>%
    filter(phase != "Burn") %>%
    group_by(hash) %>%
    summarise(
      dag = list(first(dag)), 
      score = first(score), 
      count = n()) %>%
    mutate(
      freq = count / N, 
      shdem = purrr::map_dbl(dag, ~ shd_edge_marks(true, .x)) / 2) %>%
    rename(score_0 = score)
  # Calculate scores using CS prior and ISS
  for (iss in prior_isss)
  {
    if (verbose) cat(sprintf("iss = %d\n", iss))
    column_name <- paste0("score_", iss)
    scores <- lapply(
      posterior$dag, 
      function(x) {
        score(x, df, type = "bde", prior = "cs", beta = prior, beta_iss = iss)
      })
    posterior <- posterior %>%
      mutate(!!column_name := unlist(scores))
  }
  
  # Pivot data frame to get scores in one column
  posterior <- posterior %>%
    pivot_longer(
      cols = matches("^score"),  # Select columns that start with "score"
      names_to = "iss",          # Create a new column called "iss"
      values_to = "score"  # Store the scores in a new column
    ) %>%
    mutate(
      iss = as.integer(gsub("^score_", "", iss)),
      prior = case_when(
        iss == 0 ~ "Uniform",
        TRUE ~ "CS")) %>%
    group_by(iss) %>%
    mutate(
      rank = dense_rank(desc(score)),
      normalized_score = (score - min(score)) / (max(score) - min(score))) %>%
    ungroup()
  
  posterior <- 
    posterior %>%
      mutate(experiment = experiment)
}
t_end <- Sys.time()
cat(paste0("Time elapsed: ", format_elapsed_time(t_start, t_end), "\n"))

# Closing cluster
stopCluster(cl)

# Generate a timestamp for the filename
timestamp <- format(Sys.time(), format = "%Y-%m-%d_%H-%M-%S")
# Create the filename
filename <- paste("simulations/iss/iss_", n_experiments, "_", N, "_", sample_size, "_", timestamp, ".rds", sep = "")
# Save the object as an .rds file with the timestamped filename
saveRDS(result, file = filename)
if(sound) beep(sound = "fanfare")
