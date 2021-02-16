x <- with(
  bimac,
  ouchtree(nodes=node,times=time/max(time),ancestors=ancestor,labels=species)
)

r <- paint(x,subtree=c("1"="medium","9"="large","2"="small"),
  branch=c("38"="large","2"="medium"))
plot(x,regimes=r,node.names=TRUE)

## compare to bimac['OU.LP']
h5 <- hansen(data=log(bimac['size']),tree=x,regimes=bimac['OU.LP'],
  sqrt.alpha=1,sigma=1,reltol=1e-5)
r <- paint(h5,branch=c("18"="large"),subtree=c("9"="small"))
plot(x,regimes=r,node.names=TRUE)
