library(ouch)

oubranch <- function (x0, t0, t1, alpha, sigma, theta, nstep = 10) {
  x <- x0
  t <- t0
  sigma <- t(chol(sigma))
  dt <- (t1-t0)/nstep
  for (k in 1:nstep) {
    x <- x+alpha%*%(theta-x)*dt+sigma%*%rnorm(n=ncol(sigma),sd=sqrt(dt))
  }
  x
}

ex <- function(x) {
  y <- t(chol(x))
  y[lower.tri(y,diag=T)]
}

ox <- function(x) {
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=T)] <- x
  y%*%t(y)
}

nnodes <- 501
nchar <- 2

a <- c(0.1,0.01,0.2)
alpha <- ox(a)
s <- c(0.1,0,0.1)
sigma <- ox(s)

theta <- list(
              ns=c(0,0),
              big=c(5,3),
              small=c(-2,-1)
              )

set.seed(2954189)

x <- matrix(nrow=nnodes,ncol=5,dimnames=list(NULL,c("node","ancestor","time","A","B")))
x[,1] <- 1:nnodes
node <- 1
time <- 0
x[1,'time'] <- time
while (node < nnodes) {
  x[node+1,'node'] <- node+1
  x[node+1,'time'] <- time+1
  x[node+1,'ancestor'] <- node
  x[node+2,'node'] <- node+2
  x[node+2,'time'] <- time+1
  x[node+2,'ancestor'] <- node
  node <- node+2
  time <- time+1
}
x <- as.data.frame(x)

x$reg <- as.factor(sample('ns',size=nnodes,replace=T))

x$A <- NA
x$B <- NA
x[1,c("A","B")] <- 0
for (node in 2:nnodes) {
  anc <- x[node,'ancestor']
  t0 <- x[anc,'time']
  t1 <- x[node,'time']
  r <- x[node,'reg']
  x[node,c("A","B")] <- oubranch(as.numeric(x[anc,c("A","B")]),t0,t1,alpha,sigma,theta[[r]])
}

tic <- Sys.time()
tree <- with(x,ouchtree(node,ancestor,time))
toc <- Sys.time()
print(toc-tic)

bfit <- brown(data=x[c("A","B")],tree)
print(summary(bfit))

tic <- Sys.time()
hfit <- hansen(data=x[c("A","B")],tree=tree,regimes=x['reg'],sqrt.alpha=a,sigma=s,fit=F)
toc <- Sys.time()
print(toc-tic)
print(summary(hfit))

x$reg <- as.factor(sample(c('big','small'),size=nnodes,replace=T))

x$A <- NA
x$B <- NA
x[1,c("A","B")] <- 0
for (node in 2:nnodes) {
  anc <- x[node,'ancestor']
  t0 <- x[anc,'time']
  t1 <- x[node,'time']
  r <- x[node,'reg']
  x[node,c("A","B")] <- oubranch(as.numeric(x[anc,c("A","B")]),t0,t1,alpha,sigma,theta[[r]])
}

tic <- Sys.time()
tree <- with(x,ouchtree(node,ancestor,time))
toc <- Sys.time()
print(toc-tic)

bfit <- brown(data=x[c("A","B")],tree)
print(summary(bfit))

bfit <- update(bfit,data=x[c("B","A")])
print(summary(bfit))

tic <- Sys.time()
hfit <- hansen(data=x[c("A","B")],tree=tree,regimes=x['reg'],sqrt.alpha=a,sigma=s,fit=F)
toc <- Sys.time()
print(toc-tic)
print(summary(hfit))
