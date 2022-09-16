#'
#' using LS to try to find the corrector matrix comparing Matern to Matern Neumann
#' i.e. r(t,s) = r_0(t,s) - \phi(t)^T A \phi(s)
#' where \phi, r, r_0 is known
#'
#'
source('matern.bc.R')

library(matrixcalc)

kappa <- 8
sigma <- 1
L <- 1
nu <- 5/2
n <- 50

phi <- function(t, kappa, sigma, L)
{
    m <- 2*nu+1
    r <- rep(0,m)
    for(i in 1:(m/2)){
        if(i==1){
            r[2*(i-1) + 1] <- matern.covariance(t, kappa=kappa, nu =nu, sigma = sigma )
            r[2*(i-1) + 2] <- matern.covariance(t-L, kappa=kappa, nu = nu, sigma = sigma )
        }else{
            r[2*(i-1) + 1] <- (-1)^(i-1)*matern.derivative(t, kappa=kappa, nu =nu, sigma = sigma, deriv = i-1 )
            r[2*(i-1) + 2] <- (-1)^(i-1)*matern.derivative(t-L, kappa=kappa, nu = nu, sigma = sigma, deriv = i-1 )
        }
    }
    return(r)
}

t <- L*runif(n,0,1)
s <- L*runif(n,0,1)
#matern.bc(t, t, kappa=kappa, nu = 3/2, sigma = sigma, L = L )
y_0 <- matern.covariance(abs(t-s),kappa=kappa, nu = nu, sigma = sigma)
y_n <- matern.bc(t, s, kappa=kappa, nu = nu, sigma = sigma, L = L ,bc = 'Neumann')
#
y <- y_n - y_0

# LS
# (y_i- \phi(t_i)^T A \phi(s_i) )^2 = (y_i- tr(\phi(s_i) \phi(t_i)^T )A  )^2
# = (y_i- vec(\phi(s_i) \phi(t_i)^T ) vec(A)  )^2
# = (y_i- vec(\phi(s_i) \phi(t_i)^T ) D vech(A)  )^2
D <- duplication.matrix(2*nu+1)

#A
X <- matrix(0, nrow= n, ncol= ncol(D))
for(i in 1:n){
    phi_t <- phi(t[i], kappa, sigma, L)
    phi_s <- phi(s[i], kappa, sigma, L)
    X[i,] <- t(vec((phi_s)%*%t(phi_t)))%*%D
}
E <- eigen(t(X)%*%X)
vhA <- (E$vectors%*%diag(1/pmax(E$values,10^-10))%*%t(E$vectors)) %*% (t(X)%*%y)

cat('residuals=',sum((y-X%*%vhA)^2),'\n')
A <- matrix(D%*%vhA, nrow= 2*nu+1,ncol=2*nu+1)
print(solve(A))
print(matern.derivative(0,kappa = kappa, nu=nu, sigma=sigma, deriv=4))
print(matern.derivative(L,kappa = kappa, nu=nu, sigma=sigma, deriv=4))
print(matern.derivative(0,kappa = kappa, nu=nu, sigma=sigma, deriv=2))
print(matern.derivative(L,kappa = kappa, nu=nu, sigma=sigma, deriv=2))
print(matern.derivative(L,kappa = kappa, nu=nu, sigma=sigma, deriv=1))
print(matern.covariance(L,kappa = kappa, nu=nu, sigma=sigma))
