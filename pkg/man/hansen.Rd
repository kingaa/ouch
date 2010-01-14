\name{hansen}
\alias{hansen}
\title{Hansen model of evolution along a phylogenetic tree}
\description{
  Fits the Ornstein-Uhlenbeck-based Hansen model to data.
  The fitting is done using \code{optim} or \code{subplex}.
}
\usage{
hansen(data, tree, regimes, sqrt.alpha, sigma,
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
  \item{sqrt.alpha, sigma}{
    These are used to initialize the optimization algorithm.
    The selection strength matrix \eqn{\alpha}{alpha} and the random drift variance-covariance matrix \eqn{\sigma^2}{sigma^2} are parameterized by their matrix square roots.
    Specifically, these initial guesses are each packed into lower-triangular matrices (column by column).
    The product of this matrix with its transpose is the \eqn{\alpha}{alpha} or \eqn{\sigma^2}{sigma^2} matrix.
    See Details, below.
  }
  \item{fit}{
    If \code{fit=TRUE}, then the likelihood will be maximized.
    If \code{fit=FALSE}, the likelihood will be evaluated at the specified values of \code{sqrt.alpha} and \code{sigma};
    the optima \code{theta} will be returned as well.
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
\details{
  The Hansen model for the evolution of a multivariate trait \eqn{X} along a lineage can be written as a stochastic differential equation (Ito diffusion) \deqn{dX=\alpha(\theta(t)-X(t))dt+\sigma dB(t),}{dX = alpha (theta(t)-X(t)) dt + sigma dB(t),} where \eqn{t} is time along the lineage, \eqn{\theta(t)}{theta(t)} is the optimum trait value, \eqn{B(t)} is a standard Wiener process (Brownian motion), and \eqn{\alpha}{alpha} and \eqn{\sigma}{sigma} are matrices quantifying, respectively, the strength of selection and random drift.
  Without loss of generality, one can assume \eqn{\sigma}{sigma} is lower-triangular.
  This is because only the infinitesimal variance-covariance matrix \eqn{\sigma^2=\sigma\sigma^T}{sigma^2 = sigma\%*\%transpose(sigma)} is identifiable, and for any admissible variance-covariance matrix, we can choose \eqn{\sigma}{sigma} to be lower-triangular.
  Moreover, if we view the basic model as describing evolution on a fitness landscape, then \eqn{\alpha}{alpha} will be symmetric and if we further restrict ourselves to the case of stabilizing selection, \eqn{\alpha}{alpha} will be positive definite as well.
  We make these assumptions and therefore can assume that the matrix \eqn{\alpha}{alpha} has a lower-triangular square root.

  The \code{hansen} code uses numerical optimization to maximize the likelihood.
  To do this, it parameterizes the \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices.
  Each matrix can be parameterized by \code{nchar*(nchar+1)} parameters, where \code{nchar} is the number of quantitative characters.
  Specifically, the parameters initialized by the \code{sqrt.alpha} argument of \code{hansen} are used to fill the nonzero entries of a lower-triangular matrix (in column-major order), which is then multiplied by its transpose to give the selection-strength matrix.
  The parameters specified in \code{sigma} fill the nonzero entries in the lower triangular \eqn{\sigma}{sigma} matrix.
  When \code{hansen} is executed, the numerical optimizer maximizes the likelihood over these parameters.
  The \code{print}, \code{show}, and \code{summary} methods for the resulting \code{hansentree} object display (among other things) the estimated \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices.
  The \code{coef} method extracts a named list containing not only these matrices (given as the \code{alpha.matrix} and \code{sigma.sq.matrix} elements) but also the MLEs returned by the optimizer (as \code{sqrt.alpha} and \code{sigma}, respectively).
  \strong{The latter elements should not be interpreted, but can be used to restart the algorithm, etc.}
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
