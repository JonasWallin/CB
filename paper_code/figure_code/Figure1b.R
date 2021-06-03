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

samp.time <- t.time <- samp.time2 <- samp.time.cov <- t.time <- matrix(0,ncol=n.sim,nrow=length(no))
proj <- inla.mesh.projector(mesh)

for(k in 1:length(no)){
  n.obs <- no[k]
  for(i in 1:n.sim){
    #Randomly select values of kappa and phi and construct the precision matrix
    kappa <- runif(1,1,2)
    phi <- runif(1,1,2)
    Q <- inla.spde2.precision(spde, theta = c(log(kappa), log(phi)))

    #Sample observations locations uniformly over the triangles and create the observations
    cat(k/length(no),i/n.sim,"\n")
    sample.loc <- tri_sample_simple(mesh,n.obs)
    obs.loc <- sample.loc$loc
    A <- sample.loc$A
    Y <- as.vector(A %*% X)

    #Sample using the standard GMRF method
    start_time <- Sys.time()
    R <- chol(Q,pivot=TRUE)
    reo <- attr(R, 'pivot')
    X[reo] <- solve(R,rnorm(mesh$n))
    w <- solve(t(R),t(A[,reo]))
    AQA = t(w)%*%w
    Ra <- chol(AQA,pivot=TRUE)
    reoQ <- attr(Ra, 'pivot')
    AQAv <- solve(Ra,solve(t(Ra),(A%*%X - Y)[reoQ]))
    v <- rep(0,n)
    v[reo] <- solve(R,solve(t(R),t(A[reoQ,reo])%*%AQAv))
    samp <- X - v

    samp.time[k,i] <- difftime(Sys.time(),start_time,units = "secs")

    #nSample using the new GMRF method
    start_time <- Sys.time()
    a <- 1:n.obs
    ac <- (n.obs+1):n
    SVD_Build <- c_basis2_cpp(A)
    T <- SVD_Build$T
    t.time[k,i] <- difftime(Sys.time(),start_time,units = "secs")

    start_time <- Sys.time()
    Tac <- T[ac,]
    bhat <- (t(SVD_Build$U)%*%Y)/SVD_Build$S
    muhat <- T%*%mu
    tmp <- rep(0,n-n.obs)

    Qhat <- T%*%Q%*%t(T)
    Qac <- Qhat[ac,ac]
    muac <- muhat[ac]
    R2 <- chol(Qac,pivot=TRUE)
    reoQac <- attr(R2, 'pivot')

    v <- Qhat[ac,a]%*% (bhat-muhat[a])
    mu2 <- rep(0,n-n.obs)
    R2t_inv_v <- backsolve(R2,v[reoQac],transpose=TRUE)
    mu2[reoQac] <- -backsolve(R2,R2t_inv_v)
    mu2 <- mu2 + muac
    mutilde <- t(T)%*%c(as.vector(bhat), mu2)
    tmp[reoQac] <- solve(R2,rnorm(n-n.obs))
    samp <- mutilde + t(Tac)%*%tmp

    samp.time2[k,i] <- t.time[k,i] + difftime(Sys.time(),start_time,units = "secs")

  }
}

#Compute the average time for the methods for each case
time1 <- rowMeans(samp.time)
time2 <- rowMeans(samp.time2)

#Plot the results
pdf(file="~/Dropbox/research/constraint_simulation/neurips/spde_sim.pdf", height=4, width=6)
plot(no,time1,type="l", lty=2, lwd=2, col = "white",
     ylim=c(0, max(max(time1),max(time2),max(time.cov))),
     ylab="Sampling time (s)",xlab="k")
polygon(c(no, rev(no)),
        c(apply(samp.time,1,min),rev(apply(samp.time,1,max))),
        col=rgb(0.1,0.1,0.1,alpha=0.5),
        border=NA)
lines(no,time1,lty=2,lwd=2)

polygon(c(no, rev(no)),
        c(apply(samp.time2,1,min),rev(apply(samp.time2,1,max))),
        col=rgb(1,0,0,alpha=0.5),
        border=NA)
lines(no,time2,col=2,lty=1,lwd=2)

legend(1,max(max(time1),max(time2)),legend=c("GMRF old method", "GMRF basis method"),
         col=c(1,2),lty=c(2,1),lwd=2)
dev.off()
