
library(rSPDE)
library(CB)
library(MASS)
library(Matrix)
source("matern.bc.R")
kappa <- 4





### Test with Jonas corrector
r <- function(h){(1+kappa*abs(h))*exp(-kappa*abs(h))}
dr <- function(h){kappa^2*h*exp(-kappa*abs(h))}
d2r <- function(h){-kappa^2*(1-kappa*abs(h))*exp(-kappa*abs(h))}

#cov
Make.B <- function(T){
    B11 <- matrix(c(-r(0),r(T),r(T),-r(0)),nrow=2)
    B12 <- -matrix(c(dr(0),dr(T),-dr(T),-dr(0)),nrow=2)
    B22 <- matrix(c(-d2r(0),-d2r(T),-d2r(T),-d2r(0)),nrow=2)
    return(list(B11 = B11, B12 = B12, B22 = B22))
}

phi1 <- function(t,T){
    return(c(r(t),r(T-t)))
}
phi2 <- function(t,T){
    return(c(-dr(t),dr(T-t)))
}
dphi1 <- function(t,T){
    return(c(dr(t),-dr(T-t)))
}
dphi2 <- function(t,T){
    return(c(-d2r(t),-d2r(T-t)))
}
phi <- function(t,T){
    return(c(phi1(t,T),phi2(t,T)))
}
dphi <- function(t,T){
    return(c(dphi1(t,T),dphi2(t,T)))
}

Matern32_nuem <- function(s,t,T=2, sigma=1,d=0){
    B <- Make.B(T)
    A <- solve(B$B11)
    M <- matrix(0,4,4)
    M[1,3] <-  matern.covariance(s-t,kappa=kappa,nu=3/2,sigma=sigma)
    M[1,4] <-  -matern.derivative(s-t,kappa=kappa,nu=3/2,sigma=sigma, deriv=1)
    M[2,3] <-  matern.derivative(s-t,kappa=kappa,nu=3/2,sigma=sigma, deriv=1)
    M[2,4] <- -matern.derivative(s-t, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2 )
    M[1,3] <- M[1,3]  - phi1(s,T)%*%A%*%phi1(t,T)
    M[1,4] <- M[1,4]  - phi1(s,T)%*%A%*%dphi1(t,T)
    M[2,3] <- M[2,3]  - dphi1(s,T)%*%A%*%phi1(t,T)
    M[2,4] <- M[2,4]  - dphi1(s,T)%*%A%*%dphi1(t,T)
    M <- M  + t(M)
    M[1,1] <- matern.covariance(0,kappa=kappa,nu=3/2,sigma=sigma) - phi1(s,T)%*%A%*%phi1(s,T)
    M[2,2] <-  -matern.derivative(0, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2) - dphi1(s,T)%*%A%*%dphi1(s,T)

    M[3,3] <- matern.covariance(0,kappa=kappa,nu=3/2,sigma=sigma) - phi1(t,T)%*%A%*%phi1(t,T)
    M[4,4] <- -matern.derivative(0, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2)- dphi1(t,T)%*%A%*%dphi1(t,T)

    return(M)
}

#'
#' Create covariance matrix of
#' [[X(t),X'(t)], [X(s),X'(s)]
#'
Matern32 <- function(s,t,kappa, sigma=1){



    M <- matrix(0,4,4)
    M[1,3] <-  matern.covariance(s-t,kappa=kappa,nu=3/2,sigma=sigma)
    M[1,4] <-  -matern.derivative(s-t,kappa=kappa,nu=3/2,sigma=sigma, deriv=1)
    M[2,3] <-  matern.derivative(s-t,kappa=kappa,nu=3/2,sigma=sigma, deriv=1)
    M[2,4] <- -matern.derivative(s-t, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2 )
    M <- M  + t(M)
    M[3,3] <- matern.covariance(0,kappa=kappa,nu=3/2,sigma=sigma)
    M[4,4] <-  -matern.derivative(0, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2)

    M[1,1] <- matern.covariance(0,kappa=kappa,nu=3/2,sigma=sigma)
    M[2,2] <- -matern.derivative(0, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2)
    return(M)
}


#' Create covarians matrix for
#' point 0,1,2
#'
x <- seq(0,1,length.out=4)
Sigma <- matrix(0, length(x)*2,length(x)*2)
Sigma_M <- Sigma
for(i in 1:length(x)){
    ind_i <- 2*(i-1)+1:2
    for(j in 1:length(x)){
        if(i != j){
            ind_j <- 2*(j-1)+1:2
            Sigma[c(ind_i,ind_j), c(ind_i,ind_j)] <-Matern32(x[i],x[j],kappa)
            Sigma_M[c(ind_i,ind_j), c(ind_i,ind_j)] <- Matern32_nuem(x[i],x[j], T=max(x))
        }
    }
}
#print(round((Sigma_M-Sigma),5))
#print(round((Sigma_M),4))
Q <- solve(Sigma)
Q_M <- solve(Sigma_M)
print(round((Q),4))
print(round((Q_M),4))
Q_i <- Q
diag(Q_i) <- round(diag(Q_M)/diag(Q),1)  * diag(Q_i)
A <- Matrix::sparseMatrix(i=c(1,2),
                          j = c(2,2*length(x)),
                          x=c(1,1),
                          dims=c(2,2*length(x)))
Sigma_N <- Sigma_M - t(A%*%Sigma_M)%*%solve( A%*%Sigma_M%*%t(A),(A%*%Sigma_M))

plot(x,matern.bc(x,x[1], kappa, nu = 3/2, sigma=sigma,K=40,L=max(x),bc ='Neumann'),type='l',ylab='cov neu')

lines(x,Sigma_N[1,seq(1,2*length(x),2)],col='red',lty=2)
