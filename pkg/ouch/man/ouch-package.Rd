\name{ouch-package}
\docType{package}
\alias{ouch-package}
\title{Ornstein-Uhlenbeck methods for comparative phylogenetic hypotheses}
\description{The \pkg{ouch} package provides facilities for Ornstein-Uhlenbeck based methods of phylogenetic comparative analysis.}
\section{Classes}{
  The basic class, \code{ouchtree}, is provided to encode a phylogenetic tree.
  Plot and print methods are provided.

  The class \code{browntree} derives from class \code{ouchtree} and encodes the results of fitting a Brownian Motion model to data.

  The class \code{hansentree} derives from class \code{ouchtree} and encodes the results of fitting a Hansen model to data.
}
\section{Detailed Documentation}{
  \describe{
    \item{ouchtree}{\code{\link{ouchtree}}}
    \item{Brownian motion methods}{
      \code{\link{brown}},
      \code{\link{browntree-class}}
    }
    \item{Ornstein-Uhlenbeck methods}{
      \code{\link{hansen}},
      \code{\link{hansentree-class}}
    }
    \item{simulate methods}{
      \code{\link{simulate-browntree}},
      \code{\link{simulate-hansentree}}
    }
    \item{plot methods}{
      \code{\link{plot-browntree}},
      \code{\link{plot-hansentree}}
    }
    \item{\code{ape2ouch}}{
      Convert a tree in \pkg{ape} format to \pkg{ouch} format.
      See \code{\link{ape2ouch}}
    }
  }
}
\references{
  Butler, M.A. and A.A. King (2004)
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.
  American Naturalist 164:683--695.
}
\author{Aaron A. King (kingaa at umich dot edu)}
\keyword{models}
