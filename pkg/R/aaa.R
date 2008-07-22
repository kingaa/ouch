setClass(
         'ouchtree',
         representation=representation(
           nnodes = 'integer',
           nodes = 'character',
           ancestors = 'character',
           nodelabels = 'character',
           times = 'numeric',
           root = 'integer',
           nterm = 'integer',
           term = 'integer',
           anc.numbers = 'integer',
           lineages = 'list',
           epochs = 'list',
           branch.times = 'matrix',
           depth = 'numeric'
           )
         )

setClass(
         'browntree',
         contains='ouchtree',
         representation=representation(
           call='call',
           nchar='integer',
           data='list',
           theta='list',
           sigma='numeric',
           loglik='numeric'
           )
         )

setClass(
         'hansentree',
         contains='ouchtree',
         representation=representation(
           call='call',
           nchar='integer',
           optim.diagn='list',
           hessian='matrix',
           data='list',
           regimes='list',
           beta='list',
           theta='list',
           sigma='numeric',
           alpha='numeric',
           loglik='numeric'
           )
         )

bootstrap <- function (object, nboot = 200, seed = NULL, ...) {
  stop("function 'bootstrap' is undefined for objects of class '",class(object),"'")
}
setGeneric('bootstrap')  

