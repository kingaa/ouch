% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R, R/ouchtree.R, R/brown.R, R/hansen.R
\name{print}
\alias{print}
\alias{print,ouchtree-method}
\alias{show,ouchtree-method}
\alias{show,browntree-method}
\alias{print,browntree-method}
\alias{print,hansentree-method}
\alias{show,hansentree-method}
\title{Print methods}
\usage{
\S4method{print}{ouchtree}(x, ...)

\S4method{show}{ouchtree}(object)

\S4method{show}{browntree}(object)

\S4method{print}{browntree}(x, ...)

\S4method{print}{hansentree}(x, ...)

\S4method{show}{hansentree}(object)
}
\arguments{
\item{...}{additional arguments, ignored.}

\item{object, x}{object to display.}
}
\value{
\code{print} displays the tree as a table, with (possibly) other information.

\code{print} displays the tree as a table, along with the coefficients of the fitted model and diagnostic information.

\code{print} displays the tree as a table, along with the coefficients of the fitted model and diagnostic information.
}
\description{
Print methods
}
\seealso{
Other methods: 
\code{\link{bootstrap}()},
\code{\link{coef}()},
\code{\link{logLik}},
\code{\link{plot}()},
\code{\link{simulate}()},
\code{\link{summary}()},
\code{\link{update}()}
}
\concept{methods}
