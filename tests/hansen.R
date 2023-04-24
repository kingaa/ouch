library(ouch)
set.seed(1277742405L)

tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))

h1 <- brown(log(bimac['size']),tree)
h2 <- hansen(log(bimac['size']),tree,bimac['OU.1'],sqrt.alpha=1,sigma=1)
h3 <- hansen(log(bimac['size']),tree,bimac['OU.3'],sqrt.alpha=1,sigma=1,method="Nelder-Mead")
h3 <- update(h3,method="subplex",reltol=1e-11,parscale=c(0.1,0.1))
h4 <- hansen(log(bimac['size']),tree,bimac['OU.4'],sqrt.alpha=1,sigma=1,method="L-BFGS-B")
h5 <- hansen(log(bimac['size']),tree,bimac['OU.LP'],sqrt.alpha=1,sigma=1,reltol=1e-5)
h5 <- update(h5,method='subplex',reltol=1e-11,parscale=c(0.1,0.1),hessian=TRUE)
simdat <- simulate(h5,nsim=10)
hsim <- update(h5,data=simdat[[1]])
bsim <- update(h1,data=simdat[[1]])

h2 <- update(h2,method="subplex",reltol=0)
h2 <- update(h2,method="subplex",maxit=2)

try(hansen(bimac))
try(hansen(tree=bimac))
try(hansen(tree=tree))
try(hansen(tree=tree,data=bimac$size[1:10]))
try(hansen(tree=tree,data=bimac$size))
try(hansen(tree=tree,data=setNames(bimac$size,bimac$node)))
try(hansen(tree=tree,data=setNames(bimac$size,bimac$node),sqrt.alpha=1))
try(hansen(tree=tree,data=setNames(bimac$size,bimac$node),sqrt.alpha=1,sigma=1))
try(hansen(tree=tree,data=setNames(bimac$size,bimac$node),sqrt.alpha=1,sigma=1,
  regimes=1:5))
try(hansen(tree=tree,data=setNames(bimac$size,bimac$node),sqrt.alpha=1,sigma=1,
  regimes=bimac$OU.LP))
try(hansen(tree=tree,data="bob",sqrt.alpha=1,sigma=1,
  regimes=setNames(bimac$OU.LP,bimac$node)))
try(hansen(tree=tree,data=setNames(c(bimac$size[-c(44,45)],NA,NA),bimac$node),
  sqrt.alpha=1,sigma=1,
  regimes=setNames(bimac$OU.LP,bimac$node)))
try(
  hansen(tree=tree,
    data=list(
      setNames(bimac$size,bimac$node),
      setNames(runif(n=45),bimac$node)
    ),
    sqrt.alpha=1,sigma=1,
    regimes=setNames(bimac$OU.LP,bimac$node)
  )
)
try(
  hansen(tree=tree,
    data=list(
      setNames(bimac$size,bimac$node),
      setNames(runif(n=45),bimac$node)
    ),
    sqrt.alpha=c(1,1,1),sigma=1,
    regimes=setNames(bimac$OU.LP,bimac$node)
  )
)
try(
  hansen(tree=tree,
    data=list(
      setNames(bimac$size,bimac$node),
      setNames(rnorm(n=45),bimac$node)
    ),
    sqrt.alpha=c(1,1,1),sigma=c(1,1,1),
    regimes=setNames(bimac$OU.LP[1:10],bimac$node[1:10])
  )
)
try(
  hansen(tree=tree,
    data=list(
      setNames(bimac$size,bimac$node),
      setNames(rnorm(n=45),bimac$node)
    ),
    sqrt.alpha=c(1,1,1),sigma=c(1,1,1),
    regimes=list(
      setNames(bimac$OU.LP[1:10],bimac$node[1:10]),
      setNames(bimac$OU.1,bimac$node)
    )
  )
)
try(
  hansen(tree=tree,
    data=list(
      setNames(bimac$size,bimac$node),
      setNames(rnorm(n=45),bimac$node)
    ),
    sqrt.alpha=c(1,1,1),sigma=c(1,1,1),
    regimes=list(
      as.character(setNames(bimac$OU.LP,bimac$node)),
      setNames(bimac$OU.1,bimac$node)
    )
  )
)
try(
  hansen(tree=tree,
    data=list(
      setNames(bimac$size,bimac$node),
      setNames(rnorm(n=45),bimac$node)
    ),
    sqrt.alpha=c(1,1,1),sigma=c(1,1,1),
    regimes=list(
      setNames(bimac$OU.LP,bimac$node),
      setNames(bimac$OU.1,bimac$node),
      setNames(bimac$OU.3,bimac$node)
    )
   )
)
hansen(tree=tree,
  data=list(
    setNames(bimac$size,bimac$node),
    setNames(rnorm(n=45),bimac$node)
  ),
  sqrt.alpha=c(1,1,1),sigma=c(1,1,1),
  regimes=setNames(bimac$OU.LP,bimac$node)
) -> h
hansen(tree=tree,
  data=list(
    A=setNames(bimac$size,bimac$node),
    B=setNames(rnorm(n=45),bimac$node)
  ),
  sqrt.alpha=c(1,1,1),sigma=c(1,1,1),
  regimes=setNames(bimac$OU.LP,bimac$node)
) -> h

names(as.data.frame(h))
