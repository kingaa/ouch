library(ouch)

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
