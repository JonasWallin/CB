#'
#' using LS to try to find the corrector matrix comparing Matern to Matern Neumann
#' i.e. r(t,s) = r_0(t,s) - \phi(t)^T A \phi(s)
#' where \phi, r, r_0 is known
#'
#'
library(rSPDE)
source('matern.bc.R')

library(matrixcalc)

n <- 500
K <- 100
set.seed(1)
sigma <- runif(1,1,10)
kappa <- runif(1,1,10)
L <- runif(1, 0.5, 10)
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
    r[4] <- matern.derivative(L-t, kappa=kappa, nu = 3/2, sigma = sigma, deriv=2 )
    return(r)
}

t <- L*runif(n,0,1)
s <- L*runif(n,0,1)
#matern.bc(t, t, kappa=kappa, nu = 3/2, sigma = sigma, L = L )
y_0 <- matern.covariance(abs(t-s),kappa=kappa, nu = 3/2, sigma = sigma)
y_n <- matern.bc(t, s, kappa=kappa, nu = 3/2, sigma = sigma, L = L ,bc = 'Neumann', K = K)
#
y <- y_n - y_0

# LS
# (y_i- \phi(t_i)^T A \phi(s_i) )^2 = (y_i- tr(\phi(s_i) \phi(t_i)^T )A  )^2
# = (y_i- vec(\phi(s_i) \phi(t_i)^T ) vec(A)  )^2
# = (y_i- vec(\phi(s_i) \phi(t_i)^T ) D vech(A)  )^2
D <- duplication.matrix(4)

#A
X <- matrix(0, nrow= n, ncol= ncol(D))
for(i in 1:n){
    phi_t <- phi(t[i], kappa, sigma, L)
    phi_s <- phi(s[i], kappa, sigma, L)
    X[i,] <- - t(vec((phi_s)%*%t(phi_t)))%*%D
}
vhA <- solve(t(X)%*%X, t(X)%*%y)
A <- matrix(D%*%vhA, nrow=4,ncol=4)
print('estimated:')
print(solve(A))
A.inv.theo <- A.inv(kappa,sigma = sigma, L)
A.theo <- solve(A.inv.theo)
print('theortical')
print(A.inv.theo)
#print(matern.derivative(0,kappa = kappa, nu=3/2, sigma=sigma, deriv=2))
#print(matern.derivative(L,kappa = kappa, nu=3/2, sigma=sigma, deriv=2))
#print(matern.derivative(L,kappa = kappa, nu=3/2, sigma=sigma, deriv=1))
#print(matern.covariance(L,kappa = kappa, nu=3/2, sigma=sigma))
print('diff=')
print(solve(A)- A.inv.theo)
vhA.the <- vech(A.theo)
cat('residuals=',sum((y-X%*%vhA)^2),'\n')
cat('theortical residuals=',sum((y-X%*%vhA.the)^2),'\n')
x <- seq(0,L, length.out = 100)
plot(x, matern.bc(0, x, kappa=kappa, nu = 3/2, sigma = sigma, L = L ,bc = 'Neumann',K = K),typ='l')
x.cov <- matern.covariance(x,kappa=kappa, nu = 3/2, sigma = sigma)
phi.o <-  phi(0, kappa, sigma, L)
Aphi.o <- A.theo%*%phi.o
for(i in 1:length(x)){
    phi_i <- phi(x[i], kappa, sigma, L)
    x.cov[i] <- x.cov[i] - t(Aphi.o)%*%phi_i
}
lines(x, x.cov,col='red',lty=2)

Phid01 <- cbind(phi.d(0,kappa,sigma, L),phi.d(1,kappa,sigma, L))
