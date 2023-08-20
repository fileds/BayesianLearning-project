library(igraph)
source("R/misc.R") # Loading functions for reverse edge possibility in MCMC step.

#' Perform Markov Chain Monte Carlo (MCMC) for DAG Structure Learning
#'
#' This function performs Markov Chain Monte Carlo (MCMC) to learn the structure
#' of a Directed Acyclic Graph (DAG) given data using a specified scoring function.
#'
#' @param data A data frame or matrix containing the data.
#' @param N The number of MCMC iterations.
#' @param score A scoring function used to evaluate the quality of the DAG.
#' @param start An optional starting DAG. If NULL, an empty DAG is initialized.
#' @param reverse Logical. If TRUE, considers reversing edges during MCMC steps.
#' @param verbose Logical. If TRUE, displays progress information.
#' @param ... Additional arguments passed to the scoring function.
#'
#' @return A tibble containing information about each MCMC iteration.
#'   \item{hash}{A character vector containing hash keys of the DAGs.}
#'   \item{dag}{A list of DAG structures at each iteration.}
#'   \item{probs}{A numeric vector of acceptance probabilities.}
#'   \item{scores}{A numeric vector of scores for each DAG.}
#'   \item{accepted}{An indicator vector (0 or 1) indicating whether a DAG was accepted at each iteration.}
#'
#' @examples
#' # Example usage
#' data <- data.frame(x = rnorm(100), y = rnorm(100))
#' score_function <- function(x, data) -AIC(x, data)
#' result <- dag_mcmc(data, N = 1000, score = score_function)
#'
#' @importFrom digest digest
#' @importFrom tibble tibble
#'
#' @export
dag_mcmc <- function(data, N, score, start = NULL, reverse = FALSE, 
                     verbose = TRUE, ...) {
  if (verbose && reverse) cat("\nConsidering reversing edges in steps.\n")
  if (!is.null(start))
    old <- start
  else
    old <- empty.graph(names(data))  # Initialize the starting DAG.
  # Calculate initial score.
  old_score <- do.call(score, list(x = old, data = data, ...))
  
  posterior <- lapply(1:N, function(x) list())  # Store DAGs at each iteration.
  probs <- rep(0, N)  # Store acceptance probabilities.
  scores <- rep(0, N)  # Store scores of DAGs.
  accepted <- rep(0, N)  # Store indicators of acceptance.
  hashes <- rep(character(0), N)  # Store hash keys.
  
  if (verbose) cat("\n")
  for (i in 1:N) {
    if (verbose && (i %% 10 == 0)) cat(sprintf("Iteration %4d", i), "\r")
    # MCMC step: Consider reversing edges based on the 'reverse' flag.
    if (reverse)
      new <- add_remove_or_reverse_edge(old)
    else
      new <- add_or_remove_edge(old)
    
    # Calculate new score
    new_score <- do.call(score, list(x = new, data = data, ...))
    
    # Calculate acceptance probability
    p <- acceptance_probability(old, new, old_score, new_score)
    
    # Accept or reject
    if (runif(1) < p) {
      old <- new
      old_score <- new_score
      accepted[i] <- 1
    }
    hashes[i] <- digest(amat(old))
    probs[i] <- p
    posterior[[i]] <- old
    scores[i] <- old_score
  }
  if (verbose) cat("\n")
  
  return(tibble(step = seq(1, N), hash = hashes, dag = posterior, score = scores, prob = probs, accepted = accepted))
}

#' Creates a mask where 0 indicates allowed transitions that will not exist in 'a' or cause cycles.
#'
#' @param dag The input directed acyclic graph (DAG).
#' @return A mask where 0 indicates allowed transitions that will not exist in 'a' or cause cycles.
create_mask <- function(dag) {
  a <- amat(dag)
  dist <- igraph::distances(
    igraph::graph_from_adjacency_matrix(a),
    mode = "in")
  mask <- ifelse(dist == Inf, 0, 1) + a
  
  return(mask)
}  

#' Randomly add or remove an edge in a DAG.
#' 
#' @param dag The input directed acyclic graph (DAG).
#' @return Updated DAG after adding or removing an edge.
add_or_remove_edge <- function(dag) {
  a <- amat(dag)
  n_edges <- sum(a)
  mask <- create_mask(dag)
  
  # Randomize action. Always add edge if there are no edges.
  d <- ifelse(n_edges > 1, runif(1), 1)
  
  # Randomly decide whether to add or remove an edge.
  if (d < 0.5) {
    # Remove edge
    possible_deletions <- which(a == 1)
    idx <- sample(possible_deletions, 1)
    a[idx] <- 0
  } else {
    # Add edge
    possible_additions <- which(mask == 0)
    idx <- ifelse(
      length(possible_additions) > 1,
      sample(possible_additions, 1),
      possible_additions)
    
    a[idx] <- 1
  }
  
  amat(dag) <- a
  
  return(dag)
}

#' Calculate the acceptance probability for a proposed DAG.
#' 
#' @param old The current DAG.
#' @param new The proposed DAG.
#' @param old_score The score of the current DAG.
#' @param new_score The score of the proposed DAG.

#' @return Acceptance probability.
acceptance_probability <- function(old, new, old_score, new_score) {
  # Calculate number of neighbours.
  new_neighbours <- sum(create_mask(new) == 0) + sum(amat(new))
  old_neighbours <- sum(create_mask(old) == 0) + sum(amat(old))
  
  # Calculate acceptance probability.
  p <- (old_neighbours / new_neighbours) * exp(new_score - old_score)
  
  return(min(c(1, p)))  # Ensure the probability is not greater than 1.
}