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

n.connected <- matrix(0,nrow=n.sim,ncol=length(no))

cat("Sample locations\n")
sample.loc <- list()
for(i in 1:n.sim){
    cat(i/n.sim,"\n")
    sample.loc[[i]] <- tri_sample_simple(mesh,no[length(no)])
}
cat("Sample locations done\n")


samp.time <- t.time <- samp.time2 <- samp.time.cov <- t.time <- matrix(0,ncol=n.sim,nrow=length(no))
proj <- inla.mesh.projector(mesh)

for(k in rev(1:length(no))){
    n.obs <- no[k]
    for(i in 1:n.sim){
        cat(k/length(no),i/n.sim,"\n")

        #Randomly select values of kappa and phi and construct the precision matrix
        kappa <- runif(1,1,2)
        phi <- runif(1,1,2)
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


        #Sample observations locations uniformly over the triangles and create the observations
        obs.loc <- sample.loc[[i]]$loc[1:n.obs,]
        A <- sample.loc[[i]]$A[1:n.obs,]
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
        SVD_Build <- c_basis2(A)
        T <- t(SVD_Build$T)
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
        R2t_inv_v <- solve(t(R2),v[reoQac])
        mu2[reoQac] <- -solve(R2,R2t_inv_v)
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
plot(no,time2,type="l", lty=6, lwd=2, col="white",
     ylim=c(min(c(min(time1),min(time2))),max(c(max(time1),max(time2)))),
     xlim = c(100,5350),
     ylab="Evaluation time (s)",xlab="k",log="y",yaxt="n")
axis(2, at = c(0.001,0.01,0.1,1,10,100), labels = c(0.001,0.01,0.1,1,10,100),las=2)

lines(no,time1,lty=2,lwd=2)
lines(no,time2,col=2,lty=1,lwd=2)
legend(1,max(max(time1),max(time2)),legend=c("GMRF old method", "GMRF basis method"),
       col=c(1,2),lty=c(2,1),lwd=2)
