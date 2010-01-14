\name{hansentree-methods}
\docType{methods}
\alias{hansentree-methods}
\alias{logLik,hansentree-method}
\alias{logLik-hansentree}
\alias{summary,hansentree-method}
\alias{summary-hansentree}
\alias{print,hansentree-method}
\alias{print-hansentree}
\alias{show,hansentree-method}
\alias{show-hansentree}
\alias{plot,hansentree-method}
\alias{plot-hansentree}
\alias{update,hansentree-method}
\alias{update-hansentree}
\alias{bootstrap}
\alias{bootstrap,hansentree-method}
\alias{bootstrap-hansentree}
\alias{simulate,hansentree-method}
\alias{simulate-hansentree}
\alias{coef,hansentree-method}
\alias{coef-hansentree}
\alias{as,hansentree-method}
\alias{coerce,hansentree,data.frame-method}
\title{Methods of the "hansentree" class}
\description{Methods of the "hansentree" class.}
\usage{
\S4method{logLik}{hansentree}(object)
\S4method{coef}{hansentree}(object, \dots)
\S4method{summary}{hansentree}(object, \dots)
\S4method{show}{hansentree}(object)
\S4method{print}{hansentree}(x, \dots)
\S4method{plot}{hansentree}(x, y, regimes, \dots)
\S4method{simulate}{hansentree}(object, nsim = 1, seed = NULL, \dots)
\S4method{update}{hansentree}(object, data, regimes, sqrt.alpha, sigma, \dots)
\S4method{bootstrap}{hansentree}(object, nboot = 200, seed = NULL, \dots)
\S4method{as}{hansentree}(object, class)
\S4method{coerce}{hansentree,data.frame}(from, to = "data.frame", strict = TRUE)
}
\arguments{
  \item{object}{The \code{hansentree} object.}
  \item{x}{the \code{hansentree} object.}
  \item{class}{
    character;
    name of the class to which \code{object} should be coerced.
  }
  \item{from, to}{
    the classes betwen which coercion should be performed.
  }
  \item{nsim}{
    The number of simulations to perform.
  }
  \item{nboot}{
    The number of boostraps to perform.
  }
  \item{seed}{The random seed to use in simulations.}
  \item{regimes, sqrt.alpha, sigma}{See \code{\link{hansen}}.}
  \item{data}{see \code{\link{hansen}}.}
  \item{y, strict}{Ignored.}
  \item{\dots}{
    Further arguments (either ignored or passed to underlying functions).
    In the case of \code{update}, these replace the corresponding arguments in the original call.
  }
}
\section{Methods}{
  \describe{
    \item{\code{plot()}}{
      plots the tree, with branches colored according to the selective regimes.
      See \link{plot-ouchtree} for more details.
    }
    \item{\code{print()}}{
      prints the tree as a table, along with the coefficients of the fitted model and diagnostic information.
    }
    \item{\code{show()}}{
      displays the fitted \code{hansentree} object.
    }
    \item{\code{summary()}}{
      displays information on the call, the fitted coefficients, and model selection statistics.
    }
    \item{coerce}{
      A \code{hansentree} object can be coerced to a data-frame via \code{as(object,"data.frame")}.
    }
    \item{\code{coef()}}{
      extracts the coefficients of the fitted model.
      This is a list with five elements:
      \describe{
	\item{\code{sqrt.alpha}:}{the coefficients that parameterize the alpha matrix.}
	\item{\code{sigma}:}{the coefficients that parameterize the sigma matrix.}
	\item{\code{theta}:}{
	  a list of the estimated optima, one per character.
	  Each element of the list is a vector containing one optimal value per regime.
	}
	\item{\code{alpha.matrix}:}{the alpha matrix itself.}
	\item{\code{sigma.sq.matrix}:}{the sigma-squared matrix itself.}
      }
    }
    \item{\code{logLik()}}{
      extracts the log likelihood of the fitted model.
    }
    \item{\code{update()}}{
      refines the model fit.
    }
    \item{\code{bootstrap()}}{
      performs a parametric bootstrap for confidence intervals.
    }
    \item{\code{simulate()}}{
      generates random deviates from the fitted model.
      \code{object} is the \code{hansentree} object, \code{nsim} is the desired number of replicates, and \code{seed} is (optionally) the random seed to use.
      \code{simulate} returns a list of data-frames, each comparable to the original data.
    }
  }
}
\author{Aaron A. King kingaa at umich dot edu}
\seealso{\code{\link{ouchtree}}, \code{\link{hansen}}}
\keyword{models}
