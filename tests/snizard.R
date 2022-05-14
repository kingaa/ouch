library(ouch)
if (require(ape)) {

  tree <- read.tree("snizard_tree.ape")
  squam.data <- read.csv("squamates.csv")
  tic <- Sys.time()
  tree.1 <- ape2ouch(tree)
  toc <- Sys.time()
  print(toc-tic)
  tree.data <- merge(
    as(tree.1,"data.frame"),
    squam.data,
    by.x="labels",
    by.y="Species",
    all.x=T
  )
### have a look at the times for this ostensibly ultrametric tree:
  ## print(tree.data$times)
### round these off to the nearest 0.001 to make it truly ultrametric:
  tree.data$times <- round(tree.data$times,3)
  ## print(tree.data$times)
### reconstruct the ouchtree object
  tree.2 <- with(tree.data,ouchtree(nodes,ancestors,times,labels))
  tree.data$OU2 <- as.character(tree.data$Burrow)
  tree.data$OU2[is.na(tree.data$OU2)] <- "ancestral"
  tree.data$OU2 <- as.factor(tree.data$OU2)
  ## plot(tree.1, regimes = tree.data["OU2"])
  tree.data$log.SVL.SE <- with(tree.data,log(SVL)-log(SE))
  h1 <- brown(data=tree.data['log.SVL.SE'],tree.1)
  h2 <- brown(data=tree.data['log.SVL.SE'],tree.2)
  tic <- Sys.time()
  h3 <- hansen(data=tree.data['log.SVL.SE'],tree.1,tree.data['OU2'],sqrt.alpha=0.5,sigma=1,fit=F)
  h4 <- hansen(data=tree.data['log.SVL.SE'],tree.2,tree.data['OU2'],sqrt.alpha=0.5,sigma=1,fit=F)
  toc <- Sys.time()
  print(toc-tic)
  model.fits <- c(h1,h2,h3,h4)
  names(model.fits) <- c("BM.1","BM.2","OU2.1","OU2.2")
  sapply(model.fits,function(x)c(unlist(coef(x)),summary(x)$aic.c))

}
