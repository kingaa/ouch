% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R, R/hansen.R, R/brown.R
\name{coef}
\alias{coef}
\alias{coef,hansentree-method}
\alias{coef,browntree-method}
\title{Model coefficients}
\usage{
\S4method{coef}{hansentree}(object, ...)

\S4method{coef}{browntree}(object, ...)
}
\arguments{
\item{object}{fitted model object.}

\item{...}{additional arguments, ignored.}
}
\value{
\code{coef} applied to a \code{hansentree} object returns a named list containing the estimated \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices(given as the \code{alpha.matrix} and \code{sigma.sq.matrix} elements, respectively) but also the MLE returned by the optimizer
(as \code{sqrt.alpha} and \code{sigma}, respectively).
\strong{The latter elements should not be interpreted, but can be used to restart the algorithm, etc.}

\code{coef} applied to a \code{browntree} object extracts a list with three elements:
\describe{
\item{\code{sigma}}{the coefficients of the sigma matrix.}
\item{\code{theta}}{a list of the estimated optima, one per character.}
\item{\code{sigma.sq.matrix}}{the sigma-squared matrix itself.}
}
}
\description{
\code{coef} extracts the parameters from a fitted model object.
}
\seealso{
Other methods for ouch trees: 
\code{\link{as_data_frame}},
\code{\link{bootstrap}()},
\code{\link{logLik}},
\code{\link{ouch-package}},
\code{\link{paint}()},
\code{\link{plot}()},
\code{\link{print}()},
\code{\link{simulate}()},
\code{\link{summary}()},
\code{\link{update}()}
}
\concept{methods for ouch trees}
