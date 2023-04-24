## Analysis of Anolis bimaculatus data
\donttest{ ## save time for CRAN
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),spcode))
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
plot(h5)

simdat <- simulate(h5,nsim=10)
hsim <- update(h5,data=simdat[[1]])
summary(hsim)
bsim <- update(h1,data=simdat[[1]])
summary(bsim)
}
