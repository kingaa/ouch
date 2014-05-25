\name{NEWS}
\title{News for Package 'ouch'}
\section{Changes in ouch version 2.8-3}{
  \itemize{
    \item 'as(tree,"data.frame")' now returns a data-frame with rownames identical to the node-names (tree@nodes).
  }
}
\section{Changes in ouch version 2.8-1}{
  \itemize{
    \item Functions that have been deprecated since 2.1-1 have been removed.
      The argument 'alpha' that was deprecated in favor of 'sqrt.alpha' in version 2.6-1 has been removed.
  }
}
\section{Changes in ouch version 2.7-2}{
  \itemize{
    \item Fixed typo in the documentation of 'hansen'.  Thanks to Roland Sookias for bringing this to my attention.
  }
}
\section{Changes in ouch version 2.6-1}{
  \itemize{
    \item In the 'hansen' function, the use of the argument 'alpha' was confusing.  
      	This vector parameterizes the "selection" portion of the OU model, which since version 2 has been assumed to be a positive-definite matrix (corresponding to stabilizing selection).
      	In particular, it contains the entries of the lower-triangular matrix square root of this positive definite matrix.
      To make the name more mnemonic, this argument has been renamed 'sqrt.alpha'.
      This is a superficial change meant to improve intelligibility: the algorithm itself has not changed.
      To prevent loss of backward compatibility, use of argument 'alpha' will continue to work, though it will generate a warning.
      In some future release, use of 'alpha' as an argument to 'hansen' will no longer generate a warning but will instead be passed through to the control list for 'optim' if method='optim' and will otherwise generate an error.
  }
}
\section{Changes in ouch version 2.1-1}{
  \itemize{
    \item Version 2 of ouch represents a major revision.  
      	ouch's capabilities have been substantially extended.
    \item support is now provided for multivariate characters.  
          Inference based on multivariate Ornstein-Uhlenbeck and multivariate Brownian motion models is now provided.
    \item nonultrametric trees are supported.
    \item the root character state need no longer be estimated.  
          Instead, the root state is assumed to be distributed according to the stationary distribution of the OU process.
    \item convenience functions for parametric bootstrap confidence intervals are provided.
    \item the new function 'paint' simplifies the specification of adaptive regimes in complex phylogenies.
    \item the new function 'ape2ouch' translates trees in the 'phylo' format of the ape package into 'ouchtree' format (thanks to Dan Ackerly).
    \item The new package is a complete rewrite, using S4 classes.  
      	The old functionality from ouch v. 1 remains available, but is now deprecated.
  }
}
