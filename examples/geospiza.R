### Darwin's finches.
\donttest{ ## Save time for CRAN
### The data were taken from package 'geiger' due to the latter being orphaned.
if (requireNamespace("ape")) {

  data(geospiza)
  plot(geospiza$phy)
  print(geospiza$dat)
  
### make an ouchtree out of the phy-format tree
  ot <- ape2ouch(geospiza$phy)

### merge data with tree info
  otd <- as(ot,"data.frame")
  otd <- merge(otd,geospiza$dat,by.x="labels",by.y="row.names",all=TRUE)
### row-names are used by 'hansen'
  rownames(otd) <- otd$nodes
  print(otd)
### this data-frame now contains the data as well as the tree geometry

### now remake the ouch tree
  ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
  plot(ot)

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
  plot(h1)

}
}
