\name{bn.strength class}
\alias{bn.strength class}
\alias{bn.strength-class}
\alias{bn.strength}
\title{The bn.strength class structure}
\description{

  The structure of an object of S3 class \code{bn.strength}.

}
\details{

  An object of class \code{bn.strength} is a data frame with the following
  columns (one row for each arc):

  \itemize{

    \item \code{from, to}: the nodes incident on the arc.
    \item \code{strength}: the strength of the arc. See
      \code{\link{arc.strength}}, \code{\link{boot.strength}},
      \code{\link{custom.strength}} and \code{\link{strength.plot}}
      for details.

  }

  and some additional attributes:

  \itemize{

    \item \code{nodes}: a vector of character strings, the labels of the nodes
      of the network(s) the strength were computed from.
    \item \code{method}: a character string, the method used to compute the
      strength coefficients. It can be equal to \code{test}, \code{score} or
      \code{bootstrap}.
    \item \code{threshold}: a numeric value, the threshold used to determine
      if a strength coefficient is significant.

  }

  An optional column called \code{direction} may also be present, giving the
  probability of the direction of an arc given its presence in the graph.

  Only the \code{plot()} method is defined for this class; therefore, it can
  be manipulated as a standard data frame.

}
\author{Marco Scutari}
\keyword{classes}
