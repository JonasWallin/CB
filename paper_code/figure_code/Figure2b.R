library(matlab)
library(optimx)
library(INLA)
library(CB)
library(rSPDE)
source("../util_R/nested.utils.R")
source("../util_R/build_mesh.R")


f1_base <- function(x, a=0.01, k = 1, alpha=1){
  x1x2 <- x[,1]*x[,2]
  return(k^(alpha)*exp(-a*x1x2) * x[,1] *
           (a * sin(x1x2*k ) - cos(x1x2*k)))
}
f2_base <- function(x, a=0.01, k = 1, alpha=1){
  x1x2 <- x[,1]*x[,2]
  return(k^(alpha)*exp(-a*x1x2) * x[,2]*
           (  cos(x1x2*k) - a*sin(x1x2*k )))
}

f1 <- function(x, a=0.01, ks=1, alpha=1){
  res <- rep(0,dim(x)[1])
  norm_const <- 0
  for(k in 1:ks){
    res <- f1_base(x, a, k, alpha)
  }
  return(res)
}
f2 <- function(x, a=0.01, ks=1, alpha=1){

  res <- rep(0,dim(x)[1])
  norm_const <- 0
  for(k in 1:ks){
    res <- f2_base(x, a, k, alpha)
  }
  return(res)
}

round2 <-function(x, digits=2){
  return(format(round(x,digits),nsmall=digits))
}

set.seed(23)

#Numbers of observations to perform the computations for
n.obs = round(seq(from=50,to=2000,length.out=20))

#Model parameters
sigma_Y <- 10^-4
kappa <- sqrt(8*3)/4
sigma <- 1
param_cov = c(log(sigma),log(kappa),log(sigma_Y))
param_spde <- c(param_cov[2], -param_cov[1], param_cov[3])
#Number of replicates to run to get an estimate of the variability
n.rep = 10


#Define the prediction locations and evalue the functions at those locations
n.p <- 20
x_loc_pred <- seq(from=0,to=4,length.out = n.p)
grid<-meshgrid(x_loc_pred)
pred.loc <- cbind(c(grid$x),c(grid$y))
m_1 <- f1(pred.loc)
m_2 <- f2(pred.loc)
input_data_spde_base <-build_mesh_divergence(c(0,4), 2, 60, soft_dirchlet=FALSE, soft_derivatieve=FALSE, alpha=4)

cov.time <- spde.time <-matrix(0, nrow = length(n.obs),ncol=n.rep)

for(i in 1:n.rep){
  for(k in 1:length(n.obs)){
    n_f <- n.p^2 + n.obs[k]
    #Simulate observations
    obs.loc <- 4*cbind(runif(n.obs[k]),runif(n.obs[k]))
    full.loc = rbind(pred.loc, obs.loc)
    Y <- c(f1(obs.loc),f2(obs.loc)) + sigma_Y * rnorm(2*n.obs[k])

    #Define various matrices with locations needed for the computations
    rxg <- obs.loc[,1]%*%matrix(rep(1,n.obs[k]),1,n.obs[k]) - matrix(rep(1,n.obs[k]),n.obs[k],1)%*%t(obs.loc[,1])
    ryg <- obs.loc[,2]%*%matrix(rep(1,n.obs[k]),1,n.obs[k]) - matrix(rep(1,n.obs[k]),n.obs[k],1)%*%t(obs.loc[,2])
    rxp <- obs.loc[,1]%*%matrix(rep(1,n.p^2),1,n.p^2) - matrix(rep(1,n.obs[k]),n.obs[k],1)%*%t(pred.loc[,1])
    ryp <- obs.loc[,2]%*%matrix(rep(1,n.p^2),1,n.p^2) - matrix(rep(1,n.obs[k]),n.obs[k],1)%*%t(pred.loc[,2])
    rx_full <- full.loc[,1]%*%matrix(rep(1,n_f),1,n_f) - matrix(rep(1,n_f),n_f,1)%*%t(full.loc[,1])
    ry_full <- full.loc[,2]%*%matrix(rep(1,n_f),1,n_f) - matrix(rep(1,n_f),n_f,1)%*%t(full.loc[,2])

    #Perform kriging prediction using covariance-based method
    start_time <- Sys.time()
    input_data <- list(rx = rxg, ry = ryg, rxp = rxp, ryp = ryp, y = Y, nu = 3)

    meanCust <- mean_cov(param_cov,input_data)
    cov.time[k,i] <- difftime(Sys.time(),start_time,units = "secs")

    #Perform kriging prediction using SPDE-based method
    start_time <- Sys.time()
    input_data_spde <- add_observations(Y, obs.loc, input_data_spde_base)
    Ax = inla.mesh.1d.A(input_data_spde$mesh,c(grid$x))
    Ay = inla.mesh.1d.A(input_data_spde$mesh,c(grid$y))
    Ax <- Ax[,2:(input_data_spde$n+1)]
    Ay <- Ay[,2:(input_data_spde$n+1)]
    Apred <- t(KhatriRao(t(Ax),t(Ay)))
    m_x <- CB::mean_xy(param_spde, input_data_spde,use_fake=FALSE)
    f_1_est <- Apred%*%m_x[1:input_data_spde$n^2]
    f_2_est <- Apred%*%m_x[input_data_spde$n^2+(1:input_data_spde$n^2)]
    spde.time[k,i] <- difftime(Sys.time(),start_time,units = "secs")
    cat("cov time = ", cov.time[k,i], "spde time = ", spde.time[k,i], "\n")
  }
}

#Plot results
pdf(file="~/Dropbox/research/constraint_simulation/neurips/divergence_time.pdf", height=4, width=6)
plot(n.obs,rowMeans(cov.time),type="l",lty=2,xlab="m",ylab="Prediction time (s)",col="white")
polygon(c(n.obs, rev(n.obs)),
        c(apply(cov.time,1,min),rev(apply(cov.time,1,max))),
        col=rgb(0.1,0.1,0.1,alpha=0.5),
        border=NA)
lines(n.obs,rowMeans(cov.time),lty=2,lwd=2)
polygon(c(n.obs, rev(n.obs)),
        c(apply(spde.time,1,min),rev(apply(spde.time,1,max))),
        col=rgb(1,0,0,alpha=0.5),
        border=NA)
lines(n.obs,rowMeans(spde.time),lty=1,col=2,lwd=2)
legend(1,35,legend=c("SPDE method", "Covariance-based method"),
       col=c(2,1),lty=c(1,2),lwd=2)
dev.off()
