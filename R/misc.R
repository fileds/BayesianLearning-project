# Additional functions for alternative methods. Beyond the scope of the project.


#' Create a mask to identify bidirectional and directed edges in a partially directed acyclic graph (PDAG).
#'
#' @param pdag The input partially directed acyclic graph (PDAG).
#' @return A mask with values to identify bidirectional and directed edges.
create_mec_mask <- function(pdag) {
  a <- amat(pdag)
  
  # Identify bidirected edges
  a_prime <- a + t(a)
  print(a_prime)
  
  # Identify directed edges
  a_double_prime <- a - t(a)
  print(a_double_prime)
  
  # Create mask
  mask <- 2 * as.integer(a_prime > 1) - a_double_prime
  print(mask)
} 

#' Randomly add, remove, or reverse an edge in a DAG.
#' 
#' @param dag The input directed acyclic graph (DAG).
#' @return Updated DAG after adding, removing, or reversing an edge.
add_remove_or_reverse_edge <- function(dag) {
  a <- amat(dag)  # Get adjacency matrix.
  n_edges <- sum(a)  # Count the number of edges.
  
  # Randomize action. Always add edge if there are no edges.
  d <- ifelse(n_edges > 1, runif(1), 1)
  
  if (d < 1/3) {
    # Remove edge
    possible_deletions <- which(a == 1)
    idx <- sample(possible_deletions, 1)
    a[idx] <- 0
  } else if (d > 1/3 && d < 2/3) {
    # Reverse edge
    acyc <- FALSE
    while (!acyc) {
      possible_reversals <- which(a == 1)
      idx <- sample(possible_reversals, 1)
      a <- reverse_edge(a, idx)
      tmp <- dag
      amat(tmp, check.cycles = FALSE) <- a
      acyc <- ifelse(acyclic(tmp), TRUE, FALSE)
    }
  } else {
    # Add edge
    mask <- create_mask(dag)  # Create mask of possible edge additions.
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

#' Reverse an edge in a given adjacency matrix 'a'.
#'
#' @param a The input adjacency matrix.
#' @param idx The index of the edge to reverse.
#' @return Adjacency matrix with the edge reversed.
reverse_edge <- function(a, idx) {
  n <- nrow(a)
  a[idx] <- 0
  i <- ((idx - 1) %% n) + 1
  j <- ((idx - 1) %/% n) + 1
  idx_rev <- (i - 1) * n + j
  a[idx_rev] <- 1
  return(a)
}

#' Calculate the entropy of a binary probability value.
#'
#' @param p A binary probability value.
#' @return Entropy value.
entropy <- function(p) ifelse(p > 0 & p < 1, -p * log2(p) - (1-p) * log2(1-p), 0)

#' Select a question based on the entropy of arcs in PDAGs.
#'
#' @param pdags A list of partially directed acyclic graphs (PDAGs).
#' @return A data frame with the selected question arcs.
select_question_entropy <- function(pdags) {
  amats <- lapply(pdags, function(x) as.matrix(amat(x)))
  n <- length(amats)
  freq <- Reduce("+", amats) / n
  h <- apply(freq, c(1, 2), entropy)
  from <- rownames(h)[(which(h == max(h)) - 1) %% nrow(h) + 1]
  to <- colnames(h)[ceiling(which(h == max(h)) / nrow(h))]
  
  return(data.frame(from = from, to = to))
}

#' Get the response for a given question in terms of probabilities.
#'
#' @param question A question about an edge.
#' @param true The true graph structure (DAG or PDAG).
#' @param probs Probabilities for "yes" and "no" responses (default: c(yes = 0.8, no = 0.2)).
#' @return A data frame with the response probabilities.
get_response <- function(question, true, probs = c(yes = 0.8, no = 0.2)) {
  wl <- merge(question, as.data.frame(arcs(true)))
  bl <- setdiff(question, as.data.frame(arcs(true)))
  
  answer <- cbind(wl, type = rep("w", nrow(wl)))
  answer <- rbind(answer, cbind(bl, type = rep("b", nrow(bl))))
  answer <- 
    answer %>%
    mutate(prob = case_when(
      type == "w" ~ probs[1],
      type == "b" ~ probs[2])) %>%
    select(c(from, to, prob))
  
  return(answer)
}

#' Score PDAGs
#'
#' This function calculates the scores of a list of PDAGs (Partially Directed
#' Acyclic Graphs) based on a specified scoring method. The scores quantify the
#' fit of each PDAG to the given data. Returns an error if there is no
#' consistent extension of a PDAG, see cextend in the bnlearn package.
#'
#' @param pdags A list of PDAGs.
#' @param data The input data as a data.frame.
#' @param score_ The scoring method to be used. Default is "bic".
#'
#' @return A numeric vector containing the scores for each PDAG.
#'
#' @examples
#' library(bnlearn)
#' data(learning.test)
#' df <- learning.test[1:200, ]
#' pdags <- pdag_bootstrap(df)
#' scores <- score_pdags(pdags, data = df, score_ = "bic")
#'
#' @importFrom bnlearn cextend  # Replace 'somePackage' with the actual package name
#' @importFrom bnlearn score  # Replace 'anotherPackage' with the actual package name
#'
#' @export
#'
score_pdags <- function(pdags, data, ...)
{
  scores <- c()
  for (pdag in pdags)
  {
    dag_score <- NA
    tryCatch({
      dag <- cextend(pdag)
      dag_score <- do.call(score, list(x = dag, data = data, ...))
    }, error = function(err) {
      warning("Warning: No consistent extension could be found, score set to NA.\n")
    }, finally = {
      scores <- append(scores, dag_score)
    })
  }
  
  return(scores)
}

#' Extract Unique PDAGs from a List
#'
#' This function takes a list of PDAGs and returns a list containing the unique
#' PDAGs present in the inputted list. The function also provides the counts of
#' each unique PDAG and a list of indices indicating the positions of each PDAG
#' in the initial list.
#'
#' @param pdags A list of PDAGs obtained from bootstrap resampling.
#'
#' @return A list with the unique PDAGs, their corresponding counts, and indices.
#'
#' @examples
#' pdags <- pdag_bootstrap(data)
#' unique_pdags <- extract_unique_pdags(pdags)
#' print(unique_pdags)
#'
identify_unique_pdags <- function(pdags, return_tibble = TRUE, DAG = FALSE)
{
  if (DAG) distance <- function(x, y) bnlearn::hamming(x, y)
  else distance <- function(x, y) bnlearn::shd(x, y)
  
  n_pdags <- length(pdags)
  unique_pdags <- list()
  counts <- c()
  processed <- c()
  for (i in 1:n_pdags)
  {
    if (!(i %in% processed) && i < n_pdags)
    {
      idx <- length(unique_pdags) + 1
      unique_pdags[[idx]] <- pdags[[i]]
      counts[idx] <- 1
      
      processed <- append(processed, i)
      for (j in (i+1):n_pdags)
      {
        if (!(j %in% processed) && distance(pdags[[i]], pdags[[j]]) == 0)
        {
          counts[idx] <- counts[idx] + 1
          processed <- append(processed, j)
        }
      }
    }
  }
  
  if (return_tibble)
    return(tibble(network = unique_pdags, count = counts))
  else
    return(list(pdags = unique_pdags, counts = counts, which_idx = which_idx))
}

#' Perform MCMC Over the Space of DAGs and Count Visited DAGs
#'
#' This function performs Markov Chain Monte Carlo (MCMC) simulations over the space of Directed Acyclic Graphs (DAGs). 
#' Instead of storing each DAG visited during the iterations, it counts the number of times each unique DAG was visited.
#'
#' @param data Data frame or matrix on which the DAG model is based.
#' @param N Number of MCMC iterations.
#' @param score A function to score the DAGs.
#' @param start Optional starting DAG. If NULL, an empty graph will be initialized.
#' @param reverse Logical indicating if edges should be reversed during the MCMC steps.
#' @param verbose Logical indicating if progress messages should be printed during the MCMC iterations.
#' @param ... Additional arguments to pass to the score function.
#'
#' @return A tibble with columns: 
#' \itemize{
#'   \item \code{dag} - The DAG structure.
#'   \item \code{count} - The number of times the DAG was visited during the simulation.
#'   \item \code{hash} - The hash representation of the DAG based on its adjacency matrix.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming suitable data, score function, and other required functions are defined.
#'   results <- dag_mcmc_with_count(data, 100, score_function)
#' }
dag_mcmc_with_count <- function(data, N, score, start = NULL, reverse = FALSE, 
                                verbose = TRUE, debugging = FALSE, ...) {
  if (verbose && reverse) cat("\nConsidering reversing edges in steps.\n")
  if (!is.null(start))
    old <- start
  else
    old <- empty.graph(names(data))  # Initialize the starting DAG.
  
  old_score <- do.call(score, list(x = old, data = data, ...))
  
  dag_counts <- list()  # A list to store unique DAGs and their counts.
  
  if (verbose) cat("\n")
  for (i in 1:N) {
    if (verbose && (i %% 10 == 0)) cat(sprintf("Iteration %4d", i), "\r")
    # MCMC step: Consider reversing edges based on the 'reverse' flag.
    if (reverse)
      new <- add_remove_or_reverse_edge(old)
    else
      new <- add_or_remove_edge(old)
    
    new_score <- do.call(score, list(x = new, data = data, ...))
    p <- acceptance_probability(old, new, old_score, new_score)
    if (runif(1) < p) {
      old <- new
      old_score <- new_score
    }
    
    # Hash the DAG's adjacency matrix.
    hash_key <- digest(amat(old))
    if (!is.null(dag_counts[[hash_key]])) {
      dag_counts[[hash_key]]$count <- dag_counts[[hash_key]]$count + 1
    } else {
      dag_counts[[hash_key]] <- list(dag = old, count = 1)
    }
  }
  if (verbose) cat("\n")
  
  # Extract DAGs and counts from the list
  dags <- lapply(dag_counts, `[[`, "dag")
  counts <- sapply(dag_counts, `[[`, "count")
  hashes <- names(dag_counts)
  
  # Convert into a tibble
  if (debugging)
    result_tibble <- tibble(dag = dags, count = counts, hash = hashes)
  else
    result_tibble <- tibble(dag = dags, count = counts)
  
  return(result_tibble)
}