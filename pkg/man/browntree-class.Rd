\name{browntree}
\docType{class}
\alias{browntree}
\alias{browntree-class}
\alias{logLik,browntree-method}
\alias{logLik-browntree}
\alias{summary,browntree-method}
\alias{summary-browntree}
\alias{print,browntree-method}
\alias{print-browntree}
\alias{show,browntree-method}
\alias{show-browntree}
\alias{plot,browntree-method}
\alias{plot-browntree}
\alias{simulate,browntree-method}
\alias{simulate-browntree}
\alias{update,browntree-method}
\alias{update-browntree}
\alias{bootstrap,browntree-method}
\alias{bootstrap-browntree}
\alias{coef,browntree-method}
\alias{coef-browntree}
\alias{coerce,browntree,data.frame-method}
\title{Fitted phylogenetic Brownian motion model}
\description{
  A fitted phylogenetic Brownian-motion model object.
}
\details{
  The function \code{brown} creates a \code{browntree} object by fitting
  a Brownian-motion model to data.
}
\section{Methods}{
  \describe{
    \item{\code{plot()}}{plots the tree.}
    \item{\code{print()}}{
      prints the tree as a table, along with the coefficients of the fitted model and diagnostic information.
    }
    \item{\code{show()}}{displays the fitted \code{browntree} object.}
    \item{\code{summary()}}{
      displays information on the call, the fitted coefficients, and model selection statistics.
    }
    \item{coerce}{A \code{browntree} object can be coerced to a data-frame
      via \code{as(object,"data.frame")}.}
    \item{\code{coef(object,\dots)}}{extracts the coefficients of the
      fitted model.  This is a list with three elements:
      \describe{
	\item{\code{sigma}:}{the coefficients of the sigma matrix.}
	\item{\code{theta}:}{a list of the estimated optima, one per character.}
	\item{\code{sigma..sq.matrix}:}{the sigma-squared matrix itself.}
      }
    }
    \item{\code{logLik(object,\dots)}}{
      extracts the log likelihood of the fitted model.
    }
    \item{\code{update(object, \dots)}}{
      refits the model.
      \code{object} is the \code{browntree} object.
      Additional arguments (in \code{\dots}) replace the corresponding arguments in the original call.
    }
    \item{\code{bootstrap(object, nboot = 200, seed = NULL, \dots)}}{
      performs a parametric bootstrap for estimation of confidence intervals.
      \code{object} is the \code{browntree} object.
      \code{nboot} is the number of bootstraps.
      \code{seed} allows one to fix the random seed (see \code{simulate} below).
      Additional arguments (in \code{\dots}) are passed to \code{update}.
    }
    \item{\code{simulate(object, nsim = 1, seed = NULL, \dots)}}{
      generates random deviates from the fitted model.
      \code{object} is the \code{browntree} object, \code{nsim} is the desired number of replicates, and \code{seed} is (optionally) the random seed to use.
      \code{simulate} returns a list of data-frames, each comparable to the original data.
    }
  }
}
\author{Aaron A. King kingaa at umich dot edu}
\seealso{\code{\link{ouchtree}}, \code{\link{brown}}}
\keyword{models}
