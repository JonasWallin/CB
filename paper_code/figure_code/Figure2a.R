library(matlab)
library(optimx)
library(INLA)
library(CB)
library(rSPDE)
library(doSNOW)
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

#Number of observations
n.obs = 50

#Number of basis functions in the SPDE-based methods to test
ns <- c(40,50,60,70,80,90)

#Number of replicates
n.rep = 50

#Model parameters
kappa_true <- sqrt(8*3)/4
sigma_true <- 1
sigma_Y <- 10^-4

#Prediciton locations
n.p <- 20
x_loc_pred <- seq(from=0,to=4,length.out = n.p)
grid<- matlab::meshgrid(x_loc_pred)
pred.loc <- cbind(c(grid$x),c(grid$y))
m_1 <- f1(pred.loc)
m_2 <- f2(pred.loc)
ytrue = c(m_1,m_2)
f_1 <- ytrue[1:n.p^2]
f_2 <- ytrue[n.p^2 + (1:n.p^2)]

errCust_all <- matrix(0, n.rep,2+  length(ns))
errCust_all2 <- matrix(0, n.rep, length(ns))

for(i in 1:n.rep){
  #Simulate observation locations and data
  obs.loc <- 4*cbind(runif(n.obs),runif(n.obs))
  Y <- c(f1(obs.loc),f2(obs.loc)) + sigma_Y * rnorm(2*n.obs)

  #Compute various matrices needed
  rxg <- obs.loc[,1]%*%matrix(rep(1,n.obs),1,n.obs) - matrix(rep(1,n.obs),n.obs,1)%*%t(obs.loc[,1])
  ryg <- obs.loc[,2]%*%matrix(rep(1,n.obs),1,n.obs) - matrix(rep(1,n.obs),n.obs,1)%*%t(obs.loc[,2])
  rxp <- obs.loc[,1]%*%matrix(rep(1,n.p^2),1,n.p^2) - matrix(rep(1,n.obs),n.obs,1)%*%t(pred.loc[,1])
  ryp <- obs.loc[,2]%*%matrix(rep(1,n.p^2),1,n.p^2) - matrix(rep(1,n.obs),n.obs,1)%*%t(pred.loc[,2])
  full.loc = rbind(pred.loc, obs.loc)
  n_f <- n.p^2 + n.obs
  rx_full <- full.loc[,1]%*%matrix(rep(1,n_f),1,n_f) - matrix(rep(1,n_f),n_f,1)%*%t(full.loc[,1])
  ry_full <- full.loc[,2]%*%matrix(rep(1,n_f),1,n_f) - matrix(rep(1,n_f),n_f,1)%*%t(full.loc[,2])

  input_data <- list(rx = rxg, ry = ryg, rxp = rxp, ryp = ryp, y = Y, nu = 3)

  #Estimate covariance-based model and predict
  neglik <- function(x){CB::loglike_cov(x, input_data)}
  param_cov = c(log(sigma_true),log(kappa_true),log(sigma_Y)) #start in true parameters
  param_cov = optim(param_cov, neglik)$par
  meanCust <- CB::mean_cov(param_cov,input_data)
  errCust_all[i,1] = sqrt(mean( (ytrue - mean(Y))^2 ))
  errCust_all[i,2] = sqrt(mean( (ytrue - meanCust   )^2 ))
  cat('i=',i,', rmse: ', round2(errCust_all[i,1],2),'(mean) ', round2(errCust_all[i,2],2), '(reg GP) \n',sep="")

  #Now estimate SPDE-baesed models and predict
  count <- 1
  param_spde <- c(param_cov[2], - param_cov[1], param_cov[3])
  for(n in ns){
    input_data_spde <-build_mesh_divergence(c(0,4), 4, n, 3, soft_dirchlet=FALSE,
                                            soft_derivatieve=FALSE, alpha=4)
    input_data_spde <- add_observations(Y, obs.loc, input_data_spde)

    neglik <- function(x){res = -CB::likelihood_y_gb(x, input_data_spde)
      return(res)
    }
    param_spde <- optim(param_spde, neglik, control = list(reltol=1e-3))$par

    m_x <- CB::mean_xy(param_spde, input_data_spde)

    Ax = INLA::inla.mesh.1d.A(input_data_spde$mesh,c(grid$x))
    Ay = INLA::inla.mesh.1d.A(input_data_spde$mesh,c(grid$y))
    Ax <- Ax[,2:(input_data_spde$n+1)]
    Ay <- Ay[,2:(input_data_spde$n+1)]
    Apred <- t(KhatriRao(t(Ax),t(Ay)))
    f_1_est <- Apred%*%m_x[1:input_data_spde$n^2]
    f_2_est <- Apred%*%m_x[input_data_spde$n^2+(1:input_data_spde$n^2)]

    errCust_all[i, 2+count] = sqrt(mean((ytrue- as.vector(rbind(f_1_est,f_2_est)) )^2))
    cat(round2(errCust_all[i,2+count],2),', ',sep="")
    count <- count + 1
  }
  cat('\n')
}
cat('Total rmse: ', round2(mean(errCust_all[,1]),2),'(mean) ', round2(mean(errCust_all[,2]),2), '(reg GP) ',sep="")
count <- 1
for(n in ns){
  cat(round2(mean(errCust_all[,2+count]),2),' (n=', n,') ',sep = "")
  count <- count + 1
}
cat('\n')


err1 <- f_1 - f_1_est
err2 <- f_2 - f_2_est

par(mfrow=c(1,3))
plot.vector.field(pred.loc,cbind(f_1,f_2),0.05,0.01,S2=FALSE, main = "Obs and field")
plot.vector.field(obs.loc,matrix(Y,n.obs,2),0.05,0.02,S2=FALSE, main = "vector field",col="red",add=TRUE)
plot.vector.field(pred.loc,as.matrix(cbind(f_1_est, f_2_est)),0.05,0.01,S2=FALSE, main = "reconstruction")
plot.vector.field(pred.loc,as.matrix(cbind(err1,err2)),0.05,0.01,S2=FALSE, main = "Error")

res <- colMeans(errCust_all)
sd_1 <- apply(errCust_all,2,sd)/sqrt(n.rep)
res2 <- colMeans(errCust_all2)
sd_2 <- apply(errCust_all2,2,sd)/sqrt(n.rep)


pdf(file="~/Dropbox/research/constraint_simulation/neurips/divergence_error.pdf", height=4, width=6)
par(mfrow=c(1,1))
plot(ns^2, res[3:(2+length(ns))],type='l', xlab='n',
     col='white',
     ylab = 'RMSE',
     ylim = c(0.6*min(res),1.02*max(res)),lwd=2)
polygon(c(ns[1]^2,ns[length(ns)]^2, ns[length(ns)]^2, ns[1]^2),
        c(res[2]+ 1.96*sd_1[2],res[2]+ 1.96*sd_1[2],res[2]- 1.96*sd_1[2],res[2]- 1.96*sd_1[2]),
        col=rgb(0.1,0.1,0.1,alpha=0.5),
        border=NA)

polygon(c(ns^2, rev(ns^2)),
        c(res[3:(2+length(ns))]+ 1.96*sd_1[3:(2+length(ns))],rev(res[3:(2+length(ns))]- 1.96*sd_1[3:(2+length(ns))])),
        col=rgb(1,0,0,alpha=0.5),
        border=NA)
lines(c(ns[1]^2,ns[length(ns)]^2),c(res[2],res[2]),lty=2,lwd=2)
lines(ns^2, res[3:(2+length(ns))], col='red',lwd=2)

polygon(c(ns[1]^2,ns[length(ns)]^2, ns[length(ns)]^2, ns[1]^2),
        c(res[1]+ 1.96*sd_1[1],res[1]+ 1.96*sd_1[1],res[1]- 1.96*sd_1[1],res[1]- 1.96*sd_1[1]),
        col=rgb(0,0,1,alpha=0.5),
        border=NA)
lines(c(ns[1]^2,ns[length(ns)]^2),c(res[1],res[1]),col='blue',lty=3,lwd=2)



legend(ns[length(ns)]^2 -3800,
       0.95*max(res),legend=c("SPDE method", "Mean","Covariance-based method" ),
       col=c("red", "blue",'black'), lty=c(1,3,2),lwd=2)
dev.off()
