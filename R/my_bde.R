library(bnlearn)

#' BDe Score Calculation for Bayesian Networks
#'
#' This function calculates the Bayesian Dirichlet equivalent (BDe) score for a 
#' given node in a Bayesian network.
#'
#' @param node The name of the target node for which the BDe score is calculated.
#' @param parents A character vector containing the names of parent nodes.
#' @param data A data frame containing the dataset.
#' @param args A list containing additional arguments.
#'   \code{iss} Equivalent sample size for BDe calculation.
#'   \code{pp} Parent prior value for Dirichlet prior parameters.
#'
#' @return The BDe score for the specified node and parent configuration.
#'
#' @details The BDe score is calculated based on the formula proposed by Heckerman et al.
#'
#' @examples
#' data <- data.frame(A = c('a', 'b', 'a', 'b'),
#'                    B = c('x', 'y', 'y', 'x'))
#' args <- list(iss = 100, pp = 0.5)
#' my_bde("A", c("B"), data, args)
#'
#' @seealso \code{\link{lgamma}}, \code{\link{table}}, \code{\link{configs}}
#'
#' @importFrom bnlearn configs
#'
#' @export
my_bde <- function(node, parents, data, args) {
  n = nrow(data)
  equivalent_sample_size <- args$iss
  parent_prior <- args$pp
  
  # Get the unique values for the node
  node_values <- unique(data[[node]])
  
  if (length(parents) == 0) {
    # Case with no parents
    counts <- table(data[[node]])
    alpha <- equivalent_sample_size / length(node_values)
    # First factor in Heckerman
    bde_score <- lgamma(equivalent_sample_size) - lgamma(equivalent_sample_size + n)
    # Second factor in Heckerman
    bde_score <- bde_score + sum(lgamma(alpha + counts) - lgamma(alpha))
    
  } else {
    # Case with parents
    counts <- table(data[[node]], configs(data[, parents, drop = FALSE]))
    parent_combinations <- unique(data[, parents, drop = FALSE])
    # Calculate the alpha values (Dirichlet prior parameters)
    # alpha is N_ijk^prime
    # Can we use a prior in this?
    if (!is.null(parent_prior)) {
      alpha <- equivalent_sample_size * parent_prior
    } else {
      alpha <- equivalent_sample_size / (nrow(parent_combinations) * length(node_values))
    }
    
    # Iterate through the parent combinations
    bde_score <- 0
    for (i in 1:nrow(parent_combinations)) {
      parent_combination <- as.numeric(parent_combinations[i,])
      
      # Get the counts that match the current parent combination
      subset_counts <- counts[, i]
      
      # Calculate the BDe score for the current combination
      # First factor in Heckerman
      bde_score <- bde_score + lgamma(alpha * length(node_values)) - lgamma(alpha * length(node_values) + sum(subset_counts))
      # Second factor in Heckerman
      bde_score <- bde_score + sum(lgamma(alpha + subset_counts) - lgamma(alpha))
      
    }
  }
  
  #bde_score <- bde_score + marginal_prior(node, parents)
  
  return(bde_score)
} 