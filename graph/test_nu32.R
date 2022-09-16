##
#' testing if taking matern nu=3/2
#' if we can reach it by conditinging.
##
library(rSPDE)
library(CB)
library(MASS)
library(Matrix)
source("matern.bc.R")
kappa <- 4

#phi
phi <- function(t, kappa, sigma,L){
    phi_t <- rep(0,2)
    phi_t[1] <- matern.covariance(t,kappa=kappa,nu=3/2,sigma=sigma)
    phi_t[2] <- matern.covariance(t-L,kappa=kappa,nu=3/2,sigma=sigma)
    return(phi_t)
}
#phi
phi.d1 <- function(t, kappa, sigma,L){
    phi_t <- rep(0,2)
    phi_t[1] <- matern.derivative(t,kappa=kappa,nu=3/2,sigma=sigma, deriv=1)
    phi_t[2] <- matern.derivative(t-L,kappa=kappa,nu=3/2,sigma=sigma, deriv=1)
    return(phi_t)
}
#corrector function
corrector <- function(s,t, kappa, sigma,L){

    phi_t <- phi(t, kappa, sigma, L)
    phi_s = phi(s, kappa, sigma, L)
    Ainv <- matrix(0, 2,2)
    D <- matrix(c(0,L,L,0),2,2)
    Ainv <- matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
    Ainv[1,2] <- Ainv[2,1]<- -Ainv[1,2]
    return(t(phi_t)%*%solve(Ainv,phi_s))
}
#'derivative
#' derivative with respect to t
corrector.dt <- function(s,t, kappa, sigma,L){

    phi_t <- phi.d1(t, kappa, sigma, L)
    phi_s = phi(s, kappa, sigma, L)
    Ainv <- matrix(0, 2,2)
    D <- matrix(c(0,L,L,0),2,2)
    Ainv <- matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
    Ainv[1,2] <- Ainv[2,1]<- -Ainv[1,2]
    return(t(phi_t)%*%solve(Ainv,phi_s))
}
corrector.ds <- function(s,t, kappa, sigma,L){

    phi_t <- phi(t, kappa, sigma, L)
    phi_s = phi.d1(s, kappa, sigma, L)
    Ainv <- matrix(0, 2,2)
    D <- matrix(c(0,L,L,0),2,2)
    Ainv <- matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
    Ainv[1,2] <- Ainv[2,1]<- -Ainv[1,2]
    return(t(phi_t)%*%solve(Ainv,phi_s))
}
corrector.dsdt <- function(s,t, kappa, sigma,L){

    phi_t <- phi.d1(t, kappa, sigma, L)
    phi_s = phi.d1(s, kappa, sigma, L)
    Ainv <- matrix(0, 2,2)
    D <- matrix(c(0,L,L,0),2,2)
    Ainv <- matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
    Ainv[1,2] <- Ainv[2,1]<- -Ainv[1,2]
    return(t(phi_t)%*%solve(Ainv,phi_s))
}
Matern32_der <- function(s,t,kappa=1, sigma=1, L=1){



    M <- matrix(0,4,4)
    M[1,3] <-  matern.covariance(abs(t-s),kappa=kappa,nu=3/2,sigma=sigma) + corrector(s,t,kappa,sigma,L)
    M[2,3] <- -matern.derivative(abs(t-s), kappa=kappa, nu = 3/2, sigma = sigma, deriv=1 ) + corrector.ds(s, t,kappa,sigma,L )
    M[1,4] <- matern.derivative(abs(t-s), kappa=kappa, nu = 3/2, sigma = sigma, deriv=1 ) + corrector.dt(s, t,kappa,sigma,L )
    M[2,4] <- -matern.derivative(abs(t-s), kappa=kappa, nu = 3/2, sigma = sigma, deriv=2 ) +  corrector.dsdt(s, t,kappa,sigma,L )
    M <- M  + t(M)
    M[3,3] <- matern.covariance(0,kappa=kappa,nu=3/2,sigma=sigma) + corrector(t,t,kappa,sigma,L)
    M[4,4] <-  -matern.derivative(0, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2)  +  corrector.dsdt(t, t,kappa,sigma,L )

    M[1,1] <- matern.covariance(0,kappa=kappa,nu=3/2,sigma=sigma) + corrector(s,s,kappa,sigma,L)
    M[2,2] <- -matern.derivative(0, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2)  +  corrector.dsdt(s, s,kappa,sigma,L )
    return(M)
}

# create two covarinces for [0,1], [0,1]
Sigma <- Matern32_der(0,1,kappa=kappa,L=1)
Q     <- solve(Sigma)
Q_independent <- kronecker(diag(2),Q)
Sigma_independt <- kronecker(diag(2),Sigma)

# create  covarinces for [0,2]
# for 0,1,2
Sigma_02 <- Matern32_der(0,2,kappa=kappa,L=2)
Sigma_01 <- Matern32_der(0,1,kappa=kappa,L=2)
Sigma_12 <- Matern32_der(1,2,kappa=kappa,L=2)
Sigma_012 <- matrix(0,3*2,3*2)
Sigma_012[c(1:2,5:6),c(1:2,5:6)] <- Sigma_02
Sigma_012[c(1:2,3:4),c(1:2,3:4)] <- Sigma_01
Sigma_012[c(3:4,5:6),c(3:4,5:6)] <- Sigma_12

A <- Matrix::sparseMatrix(i=c(1,1,2,2),j=c(3,5,4,6),x=c(1,-1,1,-1),dims=c(2,8))
#A <- Matrix::sparseMatrix(i=c(1,1),j=c(3,5),x=c(1,-1),dims=c(1,8))

Sigma_dep <- Sigma_independt - t(A%*%Sigma_independt)%*%solve( A%*%Sigma_independt%*%t(A),(A%*%Sigma_independt))
Sigma_dep <- Sigma_dep[-c(3,4),]
Sigma_dep <- Sigma_dep[,-c(3,4)]
matern.bc(seq(0,1,length.out=10),0,kappa,3/2,sigma,bc='Neumann',L=1)
matern.nu32(0,0,kappa,sigma,L=1)
x <- seq(0,1,length.out=100)

plot(x,matern.bc(x,0,kappa,3/2,sigma),type='l')
#print(round(solve(Sigma_012),4))
cat('constraint solve:')
print(round(Sigma_dep,3))
cat('true:')
print(round(Sigma_012,3))
cat('true:')
print(round(Sigma_012-Sigma_dep,3))
