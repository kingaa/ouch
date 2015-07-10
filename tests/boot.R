library(ouch)
data(bimac)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
h1 <- brown(
             data=log(bimac['size']),
             tree=tree
            )
h5 <- hansen(
             data=log(bimac['size']),
             tree=tree,
             regimes=bimac['OU.LP'],
             sqrt.alpha=1,
             sigma=1,
             reltol=1e-11,
             parscale=c(0.1,0.1),
             hessian=TRUE
             )

simdat <- simulate(h1,nsim=100,seed=92759587)
b1 <- sapply(simdat,function(x)summary(update(h1,data=x))$aic.c)
tic <- Sys.time()
b5 <- sapply(simdat,function(x)summary(update(h5,data=x))$aic.c)
toc <- Sys.time()
print(toc-tic)
cat("approximate 95% AIC.c cutoff",signif(quantile(b1-b5,0.95),digits=3),"\n")

boots.h1 <- bootstrap(h1,nboot=200,seed=92759587)
cat("bootstrap 95% confidence intervals for h1:\n")
print(t(sapply(boots.h1,quantile,probs=c(0.025,0.975))),digits=3)

boots.h5 <- bootstrap(h5,nboot=200,seed=92759587)
cat("bootstrap 95% confidence intervals for h5:\n")
print(t(sapply(boots.h5,quantile,probs=c(0.025,0.975))),digits=3)

