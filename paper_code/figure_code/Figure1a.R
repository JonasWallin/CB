require(fields)
require(INLA)
library(CB)
library(rSPDE)

set.seed(10)

#Vector with the numbers of observations to consider
#no <- c(10,seq(from=100,to=5000,by=100))
no <- round(seq(from=10,to=4000,length.out=20))

#Number of simulations per case
n.sim <- 10

#Create the mesh
x <- seq(from = 0, to = 10, length.out = 100)
mesh <- inla.mesh.create(lattice = inla.mesh.lattice(x = x, y = x),
                         extend = FALSE, refine = FALSE)
n <-  mesh$n

#Create the SPDE model and sample from it
spde <- inla.spde2.matern(mesh, alpha = 2)
mu <- rep(0,n)
Q <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
X <- inla.qsample(Q = Q)

like.time <- like.time2 <- like.time.cov <- t.time <- matrix(0,nrow=n.sim,ncol=length(no))

for(k in 1:length(no)){
  n.obs <- no[k]

  for(i in 1:n.sim){
    cat(k/length(no),i/n.sim,"\n")
    #Sample observations locations uniformly over the triangles and create the observations
    sample.loc <- tri_sample_simple(mesh,n.obs)
    obs.loc <- sample.loc$loc
    A <- sample.loc$A
    Y <- as.vector(A %*% X)

    #Randomly select values of kappa and phi to use for the likelihood evaluation
    kappa <- runif(1,1,2)
    phi <- runif(1,1,2)

    #Evaluate the likelihood using the covariance-based method
    start_time <- Sys.time()
    h <- as.matrix(dist(obs.loc))
    Sigma <- matern.covariance(h, kappa, 1, phi)
    R <- chol(Sigma)
    v <- solve(t(R),Y)
    like.cov <- -sum(log(diag(R))) - 0.5*t(v)%*%v
    like.time.cov[i,k] <- difftime(Sys.time(),start_time,units = "secs")

    #Construct the precision matrix needed for the GMRF methods
    start_time <- Sys.time()
    Q <- inla.spde2.precision(spde, theta = c(log(kappa), log(phi)))
    R <- chol(Q,pivot=TRUE)
    reo <- attr(R, 'pivot')
    Q.time <- difftime(Sys.time(),start_time,units = "secs")

    #Compute the likelihood using the standard GMRF method
    start_time <- Sys.time()
    w <- solve(t(R),t(A[,reo]))
    AQA = t(w)%*%w
    Ra <- chol(AQA,pivot=TRUE)
    reoQ <- attr(Ra, 'pivot')
    v <- solve(t(Ra),(Y-A%*%mu)[reoQ])
    like1 <- -sum(log(diag(Ra))) - 0.5*t(v)%*%v - n.obs*0.5*log(2*pi)

    like.time[i,k] <- Q.time + difftime(Sys.time(),start_time,units = "secs")

    #Compute the change of basis for the new method
    start_time <- Sys.time()
    a <- 1:n.obs
    ac <- (n.obs+1):n
    SVD_Build <- c_basis2_cpp(A)
    T <- SVD_Build$T
    t.time[i,k] <- difftime(Sys.time(),start_time,units = "secs")

    #Compute the likelihood using the new GMRF method
    start_time <- Sys.time()
    bhat <- (t(SVD_Build$U)%*%Y)/SVD_Build$S
    muhat <- T%*%mu

    Qhat <- T%*%Q%*%t(T)
    Qac <- Qhat[ac,ac]
    v0 <- bhat-muhat[a]
    v <- Qhat[ac,a]%*% v0

    R2 <- chol(Qac,pivot=TRUE)
    reoQac <- attr(R2, 'pivot')

    like2 <- -sum(log(SVD_Build$S)) + sum(log(diag(R))) - sum(log(diag(R2))) - n.obs*0.5*log(2*pi)
    like2 <- like2 -0.5*t(v0)%*%(Qhat[a,a]%*% v0)
    like2 <- like2 +0.5*t(v[reoQac])%*%backsolve(R2,backsolve(R2,v[reoQac],transpose=TRUE))
    like.time2[i,k] <- Q.time + t.time[i,k] + difftime(Sys.time(),start_time,units = "secs")
  }
}

#Compute the average time for the methods for each case
like.time.mean <- colMeans(like.time)
like.time.std <- apply(like.time,2,sd)/sqrt(n.sim)
like.time2.mean <- colMeans(like.time2)
like.time2.std <- apply(like.time2,2,sd)#/sqrt(n.sim)
like.time.cov.mean <- colMeans(like.time.cov)
like.time.cov.std <- apply(like.time.cov,2,sd)#/sqrt(n.sim)
t.mean <- colMeans(t.time)

#Plot the results

pdf(file="~/Dropbox/research/constraint_simulation/neurips/spde_plot.pdf", height=4, width=6)

plot(no,like.time.mean,type="l", lty=6, lwd=2, col="white",
     ylim=c(0, max(max(like.time.mean),max(like.time2.mean),max(like.time.cov.mean))),
     ylab="Evaluation time (s)",xlab="k")

polygon(c(no, rev(no)),
        c(apply(like.time,2,min),rev(apply(like.time,2,max))),
        col=rgb(0,0,1,alpha=0.5),
        border=NA)
lines(no,like.time.mean,type="l", lty=6, lwd=2, col="blue",
     ylim=c(0, max(max(like.time.mean),max(like.time2.mean),max(like.time.cov.mean))))

polygon(c(no, rev(no)),
        c(apply(like.time2,2,min),rev(apply(like.time2,2,max))),
        col=rgb(1,0,0,alpha=0.5),
        border=NA)
lines(no,like.time2.mean,col=2,lwd=2)

polygon(c(no, rev(no)),
        c(apply(like.time.cov,2,min),rev(apply(like.time.cov,2,max))),
        col=rgb(0.1,0.1,0.1,alpha=0.5),
        border=NA)
lines(no,like.time.cov.mean,col=1,lty=2,lwd=2)

legend(1,25,legend=c("GMRF basis method", "GMRF old method",  "Covariance-based method"),
       col=c(2,"blue",1),lty=c(1,6,2),lwd=2)
dev.off()
