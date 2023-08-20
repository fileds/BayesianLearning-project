library(bnlearn)
library(tidyverse)

#' @param data A data frame or matrix containing the data for which PDAG bootstrap is to be performed.
#' @param R The number of bootstrap samples to generate. Default is 200.
#' @param method The method to be used for estimating the PDAG. Default is pc.stable.
#' @param eq_class Logical indicating whether to return equivalent class CPDAGs. Default is FALSE.
#' @param ... Additional arguments to be passed to the PDAG estimation method.
#'
#' @return A list containing the PDAGs estimated from each resampled dataset.
#'
#' @examples
#' data <- read.csv("my_data.csv")
#' pdags <- pdag_bootstrap(data, R = 100, method = pc.stable, alpha = 0.5)
#'
#' @export
#'
pdag_bootstrap <- function(data, R=200, method=pc.stable, eq_class=FALSE, ...)
{
  N <- nrow(data)
  pdags <- lapply(1:R, function(x) list())
  for (i in 1:R)
  {
    sampled <- data[sample(1:N, N, replace = TRUE), ]
    pdag <- do.call(method, list(x = sampled, ...))
    if (eq_class)
      pdags[[i]] <- cpdag(pdag)
    else
      pdags[[i]] <- pdag
  }

  return(pdags)
} 