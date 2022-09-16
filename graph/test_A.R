A.inv <- function(kappa, sigma,L){

    Ainv <- matrix(0, 4,4)
    D <- matrix(c(0,L,L,0),2,2)
    A11_inv <- matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
    A11_inv[1,2] <- A11_inv[2,1] <-  -A11_inv[1,2]
    A12_inv <- matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,1)
    A22_inv <- matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,2)
    A22_inv[1,2] <- A22_inv[2,1] <-  -A22_inv[1,2]
    Ainv    <- cbind( rbind(A11_inv, A12_inv) , rbind(A12_inv, A22_inv) )
    return(Ainv)
}
