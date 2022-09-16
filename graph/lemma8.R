#'code for lemma 8

library(rSPDE)
source('matern.bc.R')

library(matrixcalc)
set.seed(3)
n <- 5000
K <- 100
kappa <- 2.2
sigma <- 1
L <- 1
A.inv <- function(kappa, sigma,L){

    Ainv <- matrix(0, 4,4)
    D <- matrix(c(0,L,L,0),2,2)
    A11_inv <- -matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
    A11_inv[1,2] <- A11_inv[2,1] <-  -A11_inv[1,2]
    A12_inv <- matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,1)
    A12_inv[1,2] <-  -A12_inv[1,2]
    A22_inv <- -matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,2)
    #A22_inv[1,2] <- A22_inv[2,1] <-  -A22_inv[1,2]
    Ainv    <- cbind( rbind(A11_inv, A12_inv) , rbind(t(A12_inv), A22_inv) )
    return(Ainv)
}


phi <- function(t, kappa, sigma, L)
{
    r <- rep(0,4)
    r[1] <- matern.covariance(t, kappa=kappa, nu = 3/2, sigma = sigma )
    r[2] <- matern.covariance(L-t, kappa=kappa, nu = 3/2, sigma = sigma )
    r[3] <- -matern.derivative(t, kappa=kappa, nu = 3/2, sigma = sigma , deriv=1)
    r[4] <- matern.derivative(L-t, kappa=kappa, nu = 3/2, sigma = sigma, deriv=1 )
    return(r)
}
phi.d <- function(t, kappa, sigma, L)
{
    r <- rep(0,4)
    r[1] <- matern.derivative(t, kappa=kappa, nu = 3/2, sigma = sigma ,deriv=1)
    r[2] <- -matern.derivative(L-t, kappa=kappa, nu = 3/2, sigma = sigma ,deriv=1)
    r[3] <- -matern.derivative(t, kappa=kappa, nu = 3/2, sigma = sigma , deriv=2)
    r[4] <- -matern.derivative(L-t, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2 )
    return(r)
}

A.theo <- solve(A.inv.theo)
A.inv.theo <- A.inv(kappa,sigma = sigma, L)
phid01 <- cbind(phi.d(0,kappa,sigma, L),phi.d(1,kappa,sigma, L))
dpAdp <- t(Phid01)%*%A.theo%*%Phid01
D <- matrix(c(0,L,L,0),2,2)
B11 <- -matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
B11[1,2] <- B11[2,1] <-  -B11[1,2]
B21 <- matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,1)
B21[1,2] <-  -B21[1,2]
B22 <- -matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,2)
C = B22

S =  B22 - t(B21)%*%solve(B11,B21)
Stilde = C - t(B21)%*%solve(B11,B21)
Atilde.inv = rbind(cbind(B11, B21), cbind(B21,-C))
dpAtildedp <- t(Phid01)%*%solve(Atilde.inv)%*%Phid01
