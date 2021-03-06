% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bimac.R
\docType{data}
\name{bimac}
\alias{bimac}
\title{Anolis bimaculatus lizard size data.}
\format{
A data frame with 45 observations on the following 8 variables.
\describe{
  \item{node}{Labels for the nodes.}
  \item{species}{Species names for extant species.}
  \item{size}{Body size (head length in mm) of extant species.}
  \item{ancestor}{Ancestral node.}
  \item{time}{Time of node.}
  \item{OU.1}{a factor with levels \code{ns}}
  \item{OU.3}{a factor with levels \code{small}, \code{medium}, \code{large}}
  \item{OU.4}{a factor with levels \code{small}, \code{medium}, \code{large}, \code{anc}}
  \item{OU.LP}{a factor with levels \code{small}, \code{medium}, \code{large}}
}
}
\source{
\Butler2004
}
\description{
This is the \emph{Anolis bimaculatus} dataset used in Butler & King (2004).
It is used to test a hypothesis of character displacement using an interspecific dataset of body sizes and current data on sympatry/allopatry.
The data frame has the following columns:
\code{species} which are species names,
\code{size} which is the phenotypic data,
and the variables \code{ancestor} and \code{time} which specify the topology of the phylogeny and the location of the nodes in time, respectively.
The columns \code{OU.1}, \code{OU.3}, \code{OU.4}, and \code{OU.LP} specify four hypothetical arrangements of selective regimes.
Explanations of the data follow:
\describe{
  \item{Body size.}{
    We use the phenotypic data and phylogeny of Losos (1990), which employed the head lengths (of males) as a proxy for body size.
    In this group of lizards, head length correlates very strongly with snout-to-vent length and the cube root of mass, which are standard measures of body size.
    The data are head lengths in mm; note that we use the log of this value in analyses.
  }
  \item{Tree topology}{
    The tree topology is encoded via two vectors: \code{ancestor} and \code{time}.
    Each node of the' phylogenetic tree has a corresponding row in the data frame, numbered from 1 to 45.
    The columns \code{ancestor} and \code{time} specify the phylogeny.
    The \code{ancestor} variable specifies the topology: it is a list indicating the ancestor of each node.
    The root node has ancestor 0.
    The variable \code{time} specifies the temporal location of each node, with the root node being at time 0.
  }
  \item{Specifications of selective regimes.}{
    (Columns \code{OU.1}, \code{OU.3}, \code{OU.4}, \code{OU.LP}).
    These columns are factors, the levels of which correspond to the \dQuote{paintings} of the respective adaptive regime hypotheses onto the phylogeny.
    Each selective regime is named (small, medium, large, etc.).
    Each column corresponds to a different painting of the selective regimes, and thus to a different hypothesis.
    In this example, there are 3 alternative models (see Butler & King 2004): \code{OU.4} is 4-regime model, \code{OU.3} is 3-regime model (all ancestors are medium), \code{OU.LP} is the linear parsimony model.
  }
}
}
\examples{
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
plot(tree,node.names=TRUE)

h1 <- brown(log(bimac['size']),tree)
h1
plot(h1)

h2 <- hansen(log(bimac['size']),tree,bimac['OU.1'],sqrt.alpha=1,sigma=1)
h2
plot(h2)

h3 <- hansen(log(bimac['size']),tree,bimac['OU.3'],sqrt.alpha=1,sigma=1)
h3
plot(h3)

h4 <- hansen(log(bimac['size']),tree,bimac['OU.4'],sqrt.alpha=1,sigma=1)
h4
plot(h4)

h5 <- hansen(log(bimac['size']),tree,bimac['OU.LP'],sqrt.alpha=1,sigma=1,reltol=1e-5)
h5 <- update(h5,method='subplex',reltol=1e-11,parscale=c(0.1,0.1),hessian=TRUE)
h5

simdat <- simulate(h5,nsim=10)
hsim <- update(h5,data=simdat[[1]])
summary(hsim)
bsim <- update(h1,data=simdat[[1]])
summary(bsim)
}
\references{
\Lazell1972

\Losos1990
}
\seealso{
Other examples: 
\code{\link{anolis.ssd}}
}
\author{
Marguerite A. Butler and Aaron A. King
}
\concept{examples}
\keyword{models}
