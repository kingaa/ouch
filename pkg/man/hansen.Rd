\name{hansen}
\alias{hansen}
\title{Hansen model of evolution along a phylogenetic tree}
\description{
  Fits the Ornstein-Uhlenbeck-based Hansen model to data.
  The fitting is done using \code{optim}.
}
\usage{
hansen(data, tree, regimes, alpha, sigma,
       fit = TRUE,
       method = c("Nelder-Mead","subplex","BFGS","L-BFGS-B"),
       hessian = FALSE, \dots)
}
\arguments{
  \item{data}{
    Phenotypic data for extant species, i.e., at the terminal twigs of the phylogenetic tree.
    This can either be a single numeric vector or a list (or data-frame) containing \code{nchar} vectors (columns) of data.
  }
  \item{tree}{
    A phylogenetic tree, specified as an \code{ouchtree} object.
  }
  \item{regimes}{
    A vector of codes, one for each node in the tree, specifying the selective regimes hypothesized to have been operative.
    Corresponding to each node, enter the code of the regime hypothesized for the branch segment terminating in that node.
    For the root node, because it has no branch segment terminating on it, the regime specification is irrelevant.
    If there are \code{nchar} quantitative characters (i.e., if \code{data} is a list with \code{nchar} columns), then \code{regimes} should be a list length 1 or \code{nchar}.
  }
  \item{alpha}{
    An initial guess of the selection strength parameters \eqn{alpha}.
    This is used to initialize the optimization algorithm.
    This vector will be translated into a symmetric matrix by the function given in \code{alpha.fn}.
  }
  \item{sigma}{
    An initial guess of the random drift parameters \eqn{sigma}, used to initialize the optimization algorithm.
    This vector is translated into a lower-triangular matrix by the function \code{sigma.fn}.
  }
  \item{fit}{
    If \code{fit=TRUE}, then the likelihood will be maximized.
    If \code{fit=FALSE}, the likelihood will be evaluated at the specified values of \code{alpha} and \code{sigma}; the optima \code{theta} will be returned as well.
  }
  \item{method}{
    The method to be used by the optimization algorithm, \code{optim}.
    See \code{\link[subplex]{subplex}} and \code{\link{optim}} for information on the available options.
  } 
  \item{hessian}{
    If \code{hessian=TRUE}, then the Hessian matrix will be computed by \code{optim}.
  }
  \item{\dots}{
    Additional arguments will be passed as \code{control} options to \code{optim} or \code{subplex}.
    See \code{\link{optim}} and \code{\link[subplex]{subplex}} for information on the available options.
  }
}
\value{
  \code{hansen} returns an object of class \code{hansentree}.
  For details on the methods of that class, see \code{\link{hansentree}}.
}
\seealso{\code{\link{ouchtree}}, \code{\link{hansentree}}, \code{\link{optim}}
}
\references{
  Butler, M.A. and A.A. King (2004)
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.
  American Naturalist 164:683-695.
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
