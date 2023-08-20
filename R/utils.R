library(digest)

# Pallette for plotting
tab10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC949", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AC")

#' Count Unique Adjacency Matrices from a List of DAGs
#'
#' This function takes a list of DAG (Directed Acyclic Graph) objects, converts each to its corresponding 
#' adjacency matrix, and then counts the number of unique adjacency matrices.
#' It returns a tibble with two columns: `dag` (containing a representative DAG 
#' for each unique adjacency matrix) and `count` (the number of times that unique matrix appears).
#'
#' @param dags A list of bn objects.
#'
#' @return A tibble with two columns: `dag` and `count`.
#' @export
#'
count_unique_dags <- function(dags) {
  
  # Convert each DAG in the list to its adjacency matrix
  matrices <- lapply(dags, amat)
  
  # Convert each matrix to its hash
  hashes <- sapply(matrices, function(m) digest(m))
  
  # Create a tibble with DAG and its hash
  dag_hash_tibble <- tibble(dag = dags, hash = hashes)
  
  # Group by hash and count each unique hash (i.e., unique DAG)
  result_tibble <- dag_hash_tibble %>%
    group_by(hash) %>%
    summarise(dag = list(first(dag)), count = n()) %>%
    ungroup()
  
  return(result_tibble)
}

#' Calculate Structural Hamming Distance (SHD) for Edge Marks
#'
#' This function calculates the Structural Hamming Distance (SHD) for edge marks
#' between a true Bayesian network structure and a discovered Bayesian network structure.
#'
#' @param true A true Bayesian network structure.
#' @param dag A discovered Bayesian network structure (DAG) to be compared against the true structure.
#'
#' @return The calculated SHD value, which represents the sum of false positives (fp) and false negatives (fn)
#'         in the comparison of skeleton and directed edge marks.
shd_edge_marks <- function(true, dag)
{
  # Compare skeletons
  skeleton_comparison <- mybnlearn::compare(mybnlearn::skeleton(true), mybnlearn::skeleton(dag))
  
  # Compare directed edges
  directed_comparison <- mybnlearn::compare(true, dag)
  
  # Calculate SHD as the sum of false positives and false negatives
  return(skeleton_comparison$fp + skeleton_comparison$fn 
         + directed_comparison$fp + directed_comparison$fn)
} 

my_shd <- function(true, dag)
{
  # Compare directed edges
  directed_comparison <- mybnlearn::compare(true, dag)
  
  # Calculate SHD as the sum of false positives and false negatives
  return(directed_comparison$fp + directed_comparison$fn)
} 

#' Generate prior probabilities from a directed acyclic graph (DAG).
#'
#' This function calculates prior probabilities for edges in a directed acyclic graph (DAG). It assigns positive and negative probabilities to edges based on user-defined probabilities.
#'
#' @param dag A bn object.
#' @param prob A numeric vector of length 2, indicating the probabilities for positive and negative edges, respectively. The default is c(0.5, 0.1).
#'
#' @return A data frame containing information about edges and their assigned probabilities.
#'
#' @import bnlearn
#' @import dplyr
#'
#' @examples
#' dag <- model2network("[A][B|A][C|B}")
#' prior <- prior_from_dag(dag, prob = c(0.7, 0.2))
#'
#' @export
prior_from_dag <- function(dag, prob = c(0.5, 0.1))
{
  # List of nodes
  nodes <- nodes(dag)
  
  # Edges with positive probability
  positives <- as.data.frame(arcs(dag))
  
  # Generate all possible combinations of edges
  negatives <- expand.grid(from = nodes, to = nodes)
  
  # Filter out self-loops (edges where "from" and "to" are the same)
  negatives <- negatives[negatives$from != negatives$to, ]
  
  # Edges with negative probability
  negatives <- anti_join(negatives, positives, by = join_by(from == from, to == to))
  
  # Assign probabilities
  negatives <- cbind(negatives, data.frame(prob = rep(prob[2], nrow(negatives))))
  positives <- cbind(positives, data.frame(prob = rep(prob[1], nrow(positives))))
  
  return(rbind(positives, negatives))
}

create_starting_dag <- function(true, n_edges) {
  start <- true
  true_arcs <- data.frame(arcs(true))
  selected_arcs <- true_arcs[sample(1:nrow(true_arcs), n_edges), ]
  arcs(start) <- selected_arcs
  
  return(start)
}

#' Create a modified graph with random edges
#'
#' This function takes an initial Bayesian network graph and randomly selects a
#' specified number of edges from its arcs to create a modified graph.
#'
#' @param true The original Bayesian network graph.
#' @param n_edges The number of edges to randomly select from the original arcs.
#'
#' @return A modified graph with randomly selected edges.
#'
#' @examples
#' true <- empty.graph(c("A", "B", "C", "D"))
#' true <- set.arc(true, "A", "B")
#' true <- set.arc(true, "A", "C")
#' modified_bn <- create_starting_dag(true, n_edges = 1)
#' print(modified_bn)
#'
#' @import bnlearn
#'
#' @export
randomize_dag_edges <- function(true, n_edges) {
  start <- true
  true_arcs <- data.frame(arcs(true))
  selected_arcs <- true_arcs[sample(1:nrow(true_arcs), n_edges), ]
  arcs(start) <- selected_arcs
  
  return(start)
}


#' Format Time Elapsed
#'
#' This function takes two time points and calculates the elapsed time between
#' them in HH:MM:SS format.
#'
#' @param t_start The starting time point.
#' @param t_end The ending time point.
#' @return A character string representing the elapsed time in HH:MM:SS format.
#' @examples
#' t_start <- Sys.time()
#' Sys.sleep(5)  # Simulating some function execution time
#' t_end <- Sys.time()
#' formatted_time <- format_elapsed_time(t_start, t_end)
#' cat(paste0("Time elapsed: ", formatted_time, "\n"))
#'
#' @export
format_elapsed_time <- function(t_start, t_end) {
  time_elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))
  hours <- floor(time_elapsed / 3600)
  minutes <- floor((time_elapsed %% 3600) / 60)
  seconds <- time_elapsed %% 60
  
  formatted_time <- sprintf("%02dh %02dm %02ds", as.integer(hours), as.integer(minutes), as.integer(seconds))
  return(formatted_time)
}
