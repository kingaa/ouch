library(ouch)

ouBirth <- read.csv("bergmann.csv",stringsAsFactors=TRUE)
tree <- with(ouBirth,ouchtree(nodes,ancestors,times,labels))

ouModelFX <- hansen(ouBirth['fakeX'],tree,ouBirth['OU2c'],sqrt.alpha=1,sigma=1,fit=FALSE)
ouModelFX@sigma <- 1
ouModelFX@theta$fakeX <- c(-2,2)
ouModelFX@sqrt.alpha <- sqrt(1)
fXlo <- simulate(ouModelFX,nsim=1,seed=1066L)
ouModelFX@sqrt.alpha <- sqrt(50)
fXhi <- simulate(ouModelFX,nsim=1,seed=1929L)

ouModelFY <- hansen(ouBirth['fakeY'],tree,ouBirth['OU2c'],sqrt.alpha=1,sigma=1,fit=FALSE)
ouModelFY@sigma <- 1
ouModelFY@theta$fakeY <- c(-2,2)
ouModelFY@sqrt.alpha <- sqrt(1)
fYlo <- simulate(ouModelFY,nsim=1,seed=1066L)
ouModelFY@sqrt.alpha <- sqrt(50)
fYhi <- simulate(ouModelFY,nsim=1,seed=1929L)

ouFakeData <- data.frame(
                         xlo=fXlo$rep.1[[1]],
                         xhi=fXhi$rep.1[[1]],
                         ylo=fYlo$rep.1[[1]],
                         yhi=fYhi$rep.1[[1]]
                         )

print(with(ouFakeData,max(abs(xlo-ylo),na.rm=TRUE)))
print(with(ouFakeData,max(abs(xhi-yhi),na.rm=TRUE)))

ouModelFXY <- hansen(
                     data=ouBirth[c("fakeX","fakeY")],
                     tree=tree,
                     regimes=ouBirth["OU2c"],
                     sqrt.alpha=c(1,0,1),
                     sigma=c(1,0,1),
                     fit=FALSE
                     )
ouModelFXY@sigma <- c(2,-3,1)
ouModelFXY@theta <- list(
                         fakeX=c(-2,2),
                         fakeY=c(-8,8)
                         )
ouModelFXY@sqrt.alpha <- c(1,0,2)
coef(ouModelFXY)

fXY <- simulate(ouModelFXY,nsim=2,seed=1929L)
plot(do.call(cbind,fXY))
