library(geiger)

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
