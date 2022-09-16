require(fields)
require(INLA)
require(CB)
require(rSPDE)

set.seed(10)

alpha <- 1

#Vector with the numbers of observations to consider
no <- c(10,seq(from=100,to=5000,by=500))
#Results in the paper use no <- c(10,seq(from=100,to=5500,by=100))

#Number of simulations per case
n.sim <- 2 #results in the paper use n.sim <- 10

#Create the mesh
x <- seq(from = 0, to = 10, length.out = 100)
mesh <- inla.mesh.create(lattice = inla.mesh.lattice(x = x, y = x),
                         extend = FALSE, refine = FALSE)
n <-  mesh$n

#Create the SPDE model and sample from it
mu <- rep(0,n)
if(alpha==1 || alpha == 2){
    spde <- inla.spde2.matern(mesh, alpha = alpha)
    Q <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
} else if(alpha == 3){
    spde <- inla.spde2.matern(mesh, alpha = 1)
    K <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
    C <- inla.mesh.fem(mesh)$c0
    Q <- K%*%C%*%K%*%C%*%K
} else if(alpha == 4){
    spde <- inla.spde2.matern(mesh, alpha = 1)
    K <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
    C <- inla.mesh.fem(mesh)$c0
    Q <- K%*%C%*%K%*%C%*%K%*%C%*%K
}
Q.nnz <- nnzero(Q)
Q <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
X <- inla.qsample(Q = Q)

like.time <- like.time2 <- like.time.cov <- t.time <- matrix(0,nrow=n.sim,ncol=length(no))
n.connected <- matrix(0,nrow=n.sim,ncol=length(no))

cat("Sample locations\n")
sample.loc <- list()
for(i in 1:n.sim){
    cat(i/n.sim,"\n")
    sample.loc[[i]] <- tri_sample_simple(mesh,no[length(no)])
}
cat("Sample locations done\n")

for(k in 1:length(no)){
    n.obs <- no[k]

    for(i in 1:n.sim){
        cat(k/length(no),i/n.sim,"\n")
        #Sample observations locations uniformly over the triangles and create the observations
        obs.loc <- sample.loc[[i]]$loc[1:n.obs,]
        A <- sample.loc[[i]]$A[1:n.obs,]
        Y <- as.vector(A %*% X)

        #Randomly select values of kappa and phi to use for the likelihood evaluation
        kappa <- runif(1,1,2)
        phi <- runif(1,1,2)

        #Evaluate the likelihood using the covariance-based method
        start_time <- Sys.time()
        h <- as.matrix(dist(obs.loc))
        Sigma <- matern.covariance(h, kappa, max(alpha-1,0.1), phi)
        R <- chol(Sigma)
        v <- solve(t(R),Y)
        like.cov <- -sum(log(diag(R))) - 0.5*t(v)%*%v
        like.time.cov[i,k] <- difftime(Sys.time(),start_time,units = "secs")

        #Construct the precision matrix needed for the GMRF methods
        start_time <- Sys.time()
        if(alpha==1 || alpha == 2){
            Q <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
        } else if(alpha == 3){
            K <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))

            KC = K%*%C
            Q <- KC%*%KC%*%K
        } else if(alpha == 4){
            K <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), 0))
            KC = K%*%C
            Q <- KC%*%KC%*%KC%*%K
        }
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
        SVD_Build <- c_basis2(A)
        T <- SVD_Build$T
        t.time[i,k] <- difftime(Sys.time(),start_time,units = "secs")
        n.connected[i,k] <- max(SVD_Build$cluster.n)
        #Compute the likelihood using the new GMRF method
        start_time <- Sys.time()
        bhat <- (t(SVD_Build$U)%*%Y)/SVD_Build$S
        muhat <- T%*%mu

        Qhat <- t(T)%*%Q%*%T
        Qac <- Qhat[ac,ac]
        v0 <- bhat-muhat[a]
        v <- Qhat[ac,a]%*% v0

        R2 <- chol(Qac,pivot=TRUE)
        reoQac <- attr(R2, 'pivot')

        like2 <- -sum(log(SVD_Build$S)) + sum(log(diag(R))) - sum(log(diag(R2))) - n.obs*0.5*log(2*pi)
        like2 <- like2 -0.5*t(v0)%*%(Qhat[a,a]%*% v0)
        like2 <- like2 +0.5*t(v[reoQac])%*%solve(R2,solve(t(R2),v[reoQac]))
        like.time2[i,k] <- Q.time + t.time[i,k] + difftime(Sys.time(),start_time,units = "secs")
    }
}

#Compute the average time for the methods for each case
like.time.mean <- colMeans(like.time)
like.time.std <- apply(like.time,2,sd)/sqrt(n.sim)
like.time2.mean <- colMeans(like.time2)
like.time2.std <- apply(like.time2,2,sd)/sqrt(n.sim)
like.time.cov.mean <- colMeans(like.time.cov)
like.time.cov.std <- apply(like.time.cov,2,sd)/sqrt(n.sim)
t.mean <- colMeans(t.time)

#Plot the results
plot(no,like.time.mean,type="l", lty=6, lwd=2, col="white",
     ylim=c(min(min(t.mean),min(like.time.mean),min(like.time2.mean),min(like.time.cov.mean)),
            max(max(like.time.mean),max(like.time2.mean),max(like.time.cov.mean))),
     xlim = c(100,5350),
     ylab="Evaluation time (s)",xlab="k",log="y",yaxt="n")
axis(2, at = c(0.001,0.01,0.1,1,10,100), labels = c(0.001,0.01,0.1,1,10,100),las=2)

lines(no,like.time.mean,type="l", lty=6, lwd=2, col="blue",
      ylim=c(0, max(max(like.time.mean),max(like.time2.mean),max(like.time.cov.mean))))

lines(no,like.time2.mean,col=2,lwd=2)

lines(no,like.time.cov.mean,col=1,lty=2,lwd=2)


lines(no,t.mean,col=2,lty=4,lwd=1)

legend(2000,0.01,
       legend=c("GMRF new", "GMRF old",
                "Covariance",
                "Basis build"),
       col=c(2,"blue",1,2),lty=c(1,6,2,3),lwd=c(2,2,2,1),
       box.lty=0,ncol=2)


plot(no,apply(n.connected,2,mean),type="l", lty=6, lwd=2, col="white",
     ylim=c(0, max(apply(n.connected,2,quantile,0.95))),
     ylab="Size of largest connected set",xlab="k")
polygon(c(no, rev(no)),
        c(apply(n.connected,2,quantile,0.05),rev(apply(n.connected,2,quantile,0.95))),
        col=rgb(0.1,0.1,0.1,alpha=0.5),
        border=NA)
lines(no,apply(n.connected,2,mean),col=1,lwd=2)
