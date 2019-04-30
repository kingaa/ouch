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
    Phenotypic data for extant species, i.e., species at the terminal twigs of the phylogenetic tree.
    This can either be a single named numeric vector, a list of \code{nchar} named vectors, or a data-frame containing \code{nchar} data variables.
    There must be an entry per variable for every node in the tree; use \code{NA} to represent missing data.
    If the data are supplied as one or more named vectors, the names attributes are taken to correspond to the node names specified when the \code{ouchtree} was constructed (see \code{\link{ouchtree}}).
    If the data are supplied as a data-frame, the rownames serve that purpose.
  }
  \item{tree}{
    A phylogenetic tree, specified as an \code{ouchtree} object.
  }
  \item{regimes}{
    A vector of codes, one for each node in the tree, specifying the selective regimes hypothesized to have been operative.
    Corresponding to each node, enter the code of the regime hypothesized for the branch segment terminating in that node.
    For the root node, because it has no branch segment terminating on it, the regime specification is irrelevant.
    If there are \code{nchar} quantitative characters, then one can specify a single set of \code{regimes} for all characters or a list of \code{nchar} regime specifications, one for each character.
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

  The \code{hansen} code uses unconstrained numerical optimization to maximize the likelihood.
  To do this, it parameterizes the \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices in a special way:
  each matrix is parameterized by \code{nchar*(nchar+1)/2} parameters, where \code{nchar} is the number of quantitative characters.
  Specifically, the parameters initialized by the \code{sqrt.alpha} argument of \code{hansen} are used to fill the nonzero entries of a lower-triangular matrix (in column-major order), which is then multiplied by its transpose to give the selection-strength matrix.
  The parameters specified in \code{sigma} fill the nonzero entries in the lower triangular \eqn{\sigma}{sigma} matrix.
  When \code{hansen} is executed, the numerical optimizer maximizes the likelihood over these parameters.
  The \code{print}, \code{show}, and \code{summary} methods for the resulting \code{hansentree} object display (among other things) the estimated \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} matrices.
  The \code{coef} method extracts a named list containing not only these matrices (given as the \code{alpha.matrix} and \code{sigma.sq.matrix} elements) but also the MLEs returned by the optimizer (as \code{sqrt.alpha} and \code{sigma}, respectively).
  \strong{The latter elements should not be interpreted, but can be used to restart the algorithm, etc.}
}
\examples{
\dontrun{
if (library(geiger,logical.return=TRUE)) {

### an example data set (Darwin's finches)
data(geospiza)
str(geospiza)
sapply(geospiza,class)

### check the correspondence between data and tree tips:
print(nc <- with(geospiza,name.check(geospiza.tree,geospiza.data)))
### looks like one of the terminal twigs has no data associated
### drop that tip:
tree <- with(geospiza,drop.tip(geospiza.tree,nc$tree_not_data))
dat <- as.data.frame(geospiza$dat)


### make an ouchtree out of the phy-format tree
ot <- ape2ouch(tree)

### merge data with tree info
otd <- as(ot,"data.frame")
### in these data, it so happens that the rownames correspond to node names
### we will exploit this correspondence in the 'merge' operation:
dat$labels <- rownames(dat)
otd <- merge(otd,dat,by="labels",all=TRUE)
rownames(otd) <- otd$nodes
print(otd)
### this data-frame now contains the data as well as the tree geometry

### now remake the ouch tree
ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

b1 <- brown(tree=ot,data=otd[c("tarsusL","beakD")])
summary(b1)

### evaluate an OU model with a single, global selective regime
otd$regimes <- as.factor("global")
h1 <- hansen(
             tree=ot,
             data=otd[c("tarsusL","beakD")],
             regimes=otd["regimes"],
             sqrt.alpha=c(1,0,1),
             sigma=c(1,0,1),
             maxit=10000
             )
summary(h1)
}
}
}
\seealso{
  \code{\link{ouchtree}}, \code{\link{hansentree}}, \code{\link{optim}}, \code{\link{bimac}}, \code{\link{anolis.ssd}}
}
\references{
  Butler, M.A. and A.A. King (2004)
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.
  American Naturalist 164:683-695.
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
