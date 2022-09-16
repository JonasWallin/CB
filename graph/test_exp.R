##
#' testing if taking expontial with nuemann boundary
#' if we can reach it by conditinging.
##

library(CB)
library(MASS)
library(Matrix)
kappa <- 1.1
exp.neumann <- function(s,t,kappa=1, L=1){
    C <- 1/sinh(kappa * L)
    return(C*(cosh( kappa * (L - abs(s-t)) ) + cosh(kappa * (s+t-L))))
}
# create two covarinces for [0,1], [0,1]
Sigma <- matrix(c(exp.neumann(0,0,kappa), exp.neumann(0,1,kappa),
                   exp.neumann(1,0,kappa), exp.neumann(1,1,kappa)),nrow=2)
Q <- solve(Sigma)
Q_independent <- kronecker(diag(2),Q)
Sigma_independt <- kronecker(diag(2),Sigma)
# create the joint for [0,1,2]
Sigma_line <- matrix(0,3,3)
points <- c(0,1,2)
for(i in 1:3){
    for(j in 1:3){
        Sigma_line[i,j] <- exp.neumann(points[i], points[j], kappa, L=2)
    }
}
Q_line <- solve(Sigma_line)
A <- Matrix::sparseMatrix(i=c(1,1),j=c(2,3),x=c(1,-1),dims=c(1,4))
Sigma_dep <- Sigma_independt - t(A%*%Sigma_independt)%*%solve( A%*%Sigma_independt%*%t(A),(A%*%Sigma_independt))

SVD_Build <- c_basis2_cpp(A)

a <- 1
ac <- 2:4
Qhat <- t(SVD_Build$T)%*%Q_independent%*%(SVD_Build$T)
Qac <- Qhat[ac,ac]
Qac_a <- Qhat[ac, a, drop=F]
Qa <- Qhat[a, a, drop=F]
Q_given <- Qac
Q0 <- (SVD_Build$T[,ac])%*%(Q_given)%*%t(SVD_Build$T[,ac,])
print(ginv(as.matrix(Sigma_dep))-Q0)
