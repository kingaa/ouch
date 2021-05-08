#' Anolis bimaculatus lizard size data
#' 
#' This is the \emph{Anolis bimaculatus} dataset used in Butler & King (2004).
#' It is used to test a hypothesis of character displacement using an interspecific dataset of body sizes and current data on sympatry/allopatry.
#' The data frame has the following columns:
#' `species` which are species names,
#' `size` which is the phenotypic data,
#' and the variables `node`, `ancestor`, and `time` which specify the topology of the phylogeny and the location of the nodes in time, respectively.
#' The columns `OU.1`, `OU.3`, `OU.4`, and `OU.LP` specify four hypothetical arrangements of selective regimes.
#' Explanations of the data follow:
#' \describe{
#'   \item{Body size.}{
#'     We use the phenotypic data and phylogeny of Losos (1990), which employed the head lengths (of males) as a proxy for body size.
#'     In this group of lizards, head length correlates very strongly with snout-to-vent length and the cube root of mass, which are standard measures of body size.
#'     The data are head lengths in mm; note that we use the log of this value in analyses.
#'   }
#'   \item{Tree structure}{
#'     The phylogenetic tree is encoded via three variables:
#'     `node`, `ancestor`, and `time`.
#'     The `node` variable gives a name to each node.
#'     The `ancestor` variable names the ancestor of each node.
#'     The root node has no ancestor (i.e., \code{ancestor=NA}).
#'     The variable `time` specifies the temporal location of each node, the root node being at time 0.
#'   }
#'   \item{Specifications of selective regimes.}{
#'     (Columns `OU.1`, `OU.3`, `OU.4`, `OU.LP`).
#'     These columns are factors, the levels of which correspond to the \dQuote{paintings} of the respective adaptive regime hypotheses onto the phylogeny (see [paint()]).
#'     Each selective regime is named (small, medium, large, etc.).
#'     Each column corresponds to a different painting of the selective regimes, and thus to a different hypothesis.
#'     In this example, there are 3 alternative models (see Butler & King 2004): `OU.4` is 4-regime model, `OU.3` is 3-regime model (all ancestors are medium), `OU.LP` is the linear parsimony model.
#'   }
#'   \item{Other variables}{
#'     In addition to the above, there is a two-letter code for each taxon (`spcode`) and the name of the island on which the taxon is found (`island`).
#'   }
#' }
#' 
#' @name bimac
#' @rdname bimac
#' @docType data
#' @family examples
#' @format A data frame with 45 observations on the following 11 variables.
#' \describe{
#'   \item{node}{Labels for the nodes.}
#'   \item{spcode}{Two-letter code for each taxon.}
#'   \item{species}{Species names for extant species.}
#'   \item{island}{Name of the island on which the population is found.}
#'   \item{size}{Body size (head length in mm) of extant species.}
#'   \item{ancestor}{Ancestral node.}
#'   \item{time}{Time of node.}
#'   \item{OU.1}{a factor with levels `ns`}
#'   \item{OU.3}{a factor with levels `small`, `medium`, `large`}
#'   \item{OU.4}{a factor with levels `small`, `medium`, `large`, `anc`}
#'   \item{OU.LP}{a factor with levels `small`, `medium`, `large`}
#' }
#' @author Marguerite A. Butler and Aaron A. King
#' @references
#' \Lazell1972
#'
#' \Losos1990
#' @source \Butler2004
#' @keywords models
#' @example examples/bimac.R
#' 
NULL
