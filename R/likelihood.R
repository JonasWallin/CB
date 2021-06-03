#'
#' computes the joint likelihood of y and b
#'
#' @param param      - (k x 1) parameter for the model
#' @param input_data - (list) $Bstar
#'                            $ystar
#'                            $T
#'                            $bstar
#'                            $alpha
#'@export
likelihood_y_b <- function(param, input_data, use_fake= FALSE){
  
  if(length(param)==7){
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- param[5:6] 
    
    sigma_Y <- exp(param[7])  
  } else if(length(param)==5) {
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- c(0,0)
    sigma_Y <- exp(param[5])
  } else {
    kappa <- exp(param[1])*c(1,1)
    tau   <- exp(param[2])*c(1,1)
    mu    <- c(0,0)
    sigma_Y <- exp(param[3])
  }
  
  
  
  if(input_data$alpha==2){
    Q_1 <- kappa[1]^4*input_data$C0 + 2*kappa[1]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_1 <-  tau[1]*Q_1
    Q_2 <- kappa[2]^4*input_data$C0 + 2*kappa[2]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_2 <- tau[2]*Q_2  
  } else if(input_data$alpha==3){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==4){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==5){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  }
  
  Q <- bdiag(Q_1,Q_2)
  
  mu_v <- c(rep(mu[1], dim(Q_1)[1]),
            rep(mu[2], dim(Q_2)[1]))
  mustar <- input_data$Tmat%*%mu_v 
  
  Qstar <-input_data$Tmat%*%Q%*% Matrix::t(input_data$Tmat)
  #Qhat <-Matrix::t(input_data$Tmat)%*%Q%*% input_data$Tmat
  
  Qstar_aac  <- Qstar[input_data$ac,input_data$a]
  Qstar_ac   <- Qstar[input_data$ac,input_data$ac]
  
  v <- (input_data$bstar-mustar[input_data$a])
  Qstar_aac_v <-  Qstar_aac%*%v
  #mu_hat <- mu[ac] - solve(Qac, Qaac%*%(bstar - mu_v[a]))
  mustar_ac <- mustar[input_data$ac]
  Qstar_ac_mu_ac <- Qstar_ac%*%mustar_ac
  # pi_{Y|AX}(y|b)
  # |Q^*|^1/2 / |hat Q|^{1/2}  (|Q^*|^1/2 - cancels with pi_{AX}(b))
  Qhathat <- Qstar_ac + sigma_Y^(-2) * input_data$BtB
  if(use_fake){
    Qhathat <- Qhathat + input_data$BtB_fake
  }
  
  mu_hathat <- Qstar_ac_mu_ac -  Qstar_aac_v + sigma_Y^(-2) * input_data$BtY
   
   
  R_hathat <- chol(Qhathat,pivot=TRUE)
  reo <- attr(R_hathat, 'pivot')
  #(Qtilde%*%mu_tilde)^T Qtilde^-1  (Qtilde%*%mu_tilde)
  w <- Matrix::solve(Matrix::t(R_hathat),mu_hathat[reo])
  lik <-  0.5 * (Matrix::t(w)%*%w)
  # mu_hat^T Q_ac_ac mu_hat
  # mu_hat = mu_ac + Q_{ac_ac}^{-1} Q_{ac_a}(b-\mu_a)
  # terms cancel below
  lik <- lik - 0.5 * (t(mustar_ac)%*%Qstar_ac_mu_ac)
  lik <- lik +       (t(mustar_ac)%*%Qstar_aac_v)
  #lik <- lik - 0.5 * (Matrix::t(w)%*%w)
  
  
  # y*^T * y
  lik <- lik - sigma_Y^(-2) * 0.5 * input_data$ysTys
  # 1/sigma_Y^m
  lik <- lik - length(input_data$y) * log(sigma_Y)
  # #-0.5*log(|Q_hathat|)
  lik <- lik - sum(log(diag(R_hathat)))
  
  #pi_{AX}(b)
  #R_ac <- chol(Qac,pivot=TRUE)
  #reo <- attr(R_ac, 'pivot')
  
  # (b^* -\mu^*)^T \tilde{Q}(b^* -\mu^*)^T 
  lik <- lik -0.5 * (Matrix::t(v)%*%Qstar[input_data$a, input_data$a]%*%v)
  # canceled with terms above
  #w <- Matrix::solve(Matrix::t(R_ac), Qaacv[reo])
  #lik <- lik + 0.5 * (Matrix::t(w)%*%w)
  # lik det Q
  R_1 <- chol(Q_1,pivot=TRUE)
  R_2 <- chol(Q_2,pivot=TRUE)
  lik <- lik + sum(log(diag(R_1)))
  lik <- lik + sum(log(diag(R_2)))
  

  
  return(as.vector(lik ))
}


#'
#' Computes the loglikelihood (y,b) using the standard likelihood formulation
#' 
#' @param param      - (k x 1) parameter for the model
#' @param input_data - (list) $Bstar
#'                            $ystar
#'                            $T
#'                            $bstar
#'                            $alpha
#'@export                          
likelihood_y_b2<-function(param, input_data){
  if(length(param)==7){
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- param[5:6] 
    
    sigma_Y <- exp(param[7])  
  } else if(length(param)==5) {
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- c(0,0)
    sigma_Y <- exp(param[5])
  } else {
    kappa <- exp(param[1])*c(1,1)
    tau   <- exp(param[2])*c(1,1)
    mu    <- c(0,0)
    sigma_Y <- exp(param[3])
  }
  
  
  if(input_data$alpha==2){
    Q_1 <- kappa[1]^4*input_data$C0 + 2*kappa[1]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_1 <-  tau[1]*Q_1
    Q_2 <- kappa[2]^4*input_data$C0 + 2*kappa[2]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_2 <- tau[2]*Q_2  
  } else if(input_data$alpha==3){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==4){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==5){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  }
  Sigma_1 <- solve(Q_1)
  Sigma_2 <- solve(Q_2)
  Sigma <- bdiag(Sigma_1,Sigma_2)
  # likelihood A_conX,A_obsX
  #
  Sigma_Y <- input_data$A_obs%*%Sigma%*%Matrix::t(input_data$A_obs)
  Matrix::diag(Sigma_Y) <- Matrix::diag(Sigma_Y) + sigma_Y^2
  Sigma_b <- input_data$A_con%*%Sigma%*%Matrix::t(input_data$A_con)
  Cov_Yb  <- input_data$A_obs%*%Sigma%*%Matrix::t(input_data$A_con)
  
  Sigma_Yb <- rbind( cbind(Sigma_Y, Cov_Yb),
                     cbind(Matrix::t(Cov_Yb), Sigma_b))
  mu_v <- c(rep(mu[1], dim(Q_1)[1]),
            rep(mu[2], dim(Q_2)[1]))
  
  
  mu_Yb <- as.vector(rbind(input_data$A_obs%*%mu_v,input_data$A_con%*%mu_v))
  R <- chol(forceSymmetric(Sigma_Yb))
  w <- Matrix::solve(Matrix::t(R), c(input_data$y, input_data$b)-mu_Yb)
  lik <- -(0.5 * t(w)%*%w) - sum(log(diag(R)))
  return(as.vector(lik))
}



##
#' computes posterior mean
#'
#'
##
mean_xy <- function(param, input_data, use_fake=FALSE){
  if(length(param)==7){
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- param[5:6] 
    
    sigma_Y <- exp(param[7])  
  } else if(length(param)==5) {
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- c(0,0)
    sigma_Y <- exp(param[5])
  } else {
    kappa <- exp(param[1])*c(1,1)
    tau   <- exp(param[2])*c(1,1)
    mu    <- c(0,0)
    sigma_Y <- exp(param[3])
  }
  
  

  if(input_data$alpha==2){
    Q_1 <- kappa[1]^4*input_data$C0 + 2*kappa[1]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_1 <-  tau[1]*Q_1
    Q_2 <- kappa[2]^4*input_data$C0 + 2*kappa[2]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_2 <- tau[2]*Q_2  
  } else if(input_data$alpha==3){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==4){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==5){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  }

  
  Q <- bdiag(Q_1,Q_2)
  
  mu_v <- c(rep(mu[1], dim(Q_1)[1]),
            rep(mu[2], dim(Q_2)[1]))
  mus <- input_data$Tmat%*%mu_v 
  
  Qs <-input_data$Tmat%*%Q%*% Matrix::t(input_data$Tmat)

  Qac <- Qs[input_data$ac,input_data$ac]
  
  Qaacv <-  Qs[input_data$ac,input_data$a]%*%(input_data$bstar-mus[input_data$a])
  Qac_mu_ac <- Qac%*%mus[input_data$ac]
  
  
  mu_tilde <- Qac_mu_ac -  Qaacv + sigma_Y^(-2) * input_data$BtY
  
  Qtilde <- Qac + sigma_Y^(-2) * input_data$BtB
  if(use_fake){
    Qtilde <- Qtilde + input_data$BtB_fake
  }
  #R_tilde <- chol(Qtilde,pivot=TRUE)
  #reo <- attr(R_tilde, 'pivot')
  #mean_ac <- rep(0,length(mu_tilde)) 
  #mean_ac[reo] <- Matrix::solve(R_tilde,Matrix::solve(Matrix::t(R_tilde),mu_tilde[reo]))
  L_tilde <- Matrix::Cholesky(Qtilde,perm=TRUE,super=T,LDL = F)
  mean_ac <- Matrix::solve(L_tilde,mu_tilde,system='A')
  
  mean_X <-  as.vector(Matrix::t(input_data$Tmat)%*%c(as.vector(input_data$bstar), as.vector(mean_ac)))
  return(mean_X)
}

#'
#' computes the joint likelihood of y and b
#'
#' @param param      - (k x 1) parameter for the model
#' @param input_data - (list) $Bstar
#'                            $ystar
#'                            $T
#'                            $bstar
#'                            $alpha
#'@export
likelihood_y_gb <- function(param, input_data, use_fake= FALSE){
  
  if(length(param)==7){
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- param[5:6] 
    
    sigma_Y <- exp(param[7])  
  } else if(length(param)==5) {
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- c(0,0)
    sigma_Y <- exp(param[5])
  } else {
    kappa <- exp(param[1])*c(1,1)
    tau   <- exp(param[2])*c(1,1)
    mu    <- c(0,0)
    sigma_Y <- exp(param[3])
  }
  
  
  
  if(input_data$alpha==2){
    Q_1 <- kappa[1]^4*input_data$C0 + 2*kappa[1]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_1 <-  tau[1]*Q_1
    Q_2 <- kappa[2]^4*input_data$C0 + 2*kappa[2]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_2 <- tau[2]*Q_2  
  } else if(input_data$alpha==3){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==4){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==5){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  }
  
  Q <- bdiag(Q_1,Q_2)
  
  mu_v <- c(rep(mu[1], dim(Q_1)[1]),
            rep(mu[2], dim(Q_2)[1]))
  mustar <- input_data$Tmat%*%mu_v 
  
  Qstar <-input_data$Tmat%*%Q%*% Matrix::t(input_data$Tmat)
  #Qhat <-Matrix::t(input_data$Tmat)%*%Q%*% input_data$Tmat
  
  Qstar_aac  <- Qstar[input_data$ac,input_data$a]
  Qstar_ac   <- Qstar[input_data$ac,input_data$ac]
  
  v <- (input_data$bstar-mustar[input_data$a])
  Qstar_aac_v <-  Qstar_aac%*%v
  #mu_hat <- mu[ac] - solve(Qac, Qaac%*%(bstar - mu_v[a]))
  mustar_ac <- mustar[input_data$ac]
  Qstar_ac_mu_ac <- Qstar_ac%*%mustar_ac
  # pi_{Y|AX}(y|b)
  # |Q^*|^1/2 / |hat Q|^{1/2}  (|Q^*|^1/2 - cancels with pi_{AX}(b))
  Qhathat <- Qstar_ac + sigma_Y^(-2) * input_data$BtB
  if(use_fake){
    Qhathat <- Qhathat + input_data$BtB_fake
  }
  
  mu_hathat <- Qstar_ac_mu_ac -  Qstar_aac_v + sigma_Y^(-2) * input_data$BtY
  #(Qtilde%*%mu_tilde)^T Qtilde^-1  (Qtilde%*%mu_tilde)
  #R_ac <- chol(Qhathat, pivot=TRUE)
  #reo <- attr(R_ac, 'pivot')
  #w <- Matrix::solve(t(R_ac),mu_hathat[reo])
  L_hathat <- Matrix::Cholesky(Qhathat,perm=TRUE,super=T,LDL = F)
  w <- Matrix::solve(L_hathat,Matrix::solve(L_hathat,mu_hathat,system='P'),system='L')
  lik <-  0.5 * (Matrix::t(w)%*%w)
  # #-0.5*log(|Q_hathat|)
  lik <- lik - Matrix::determinant(L_hathat)$mod
  #lik <- lik - sum(log(diag(R_ac)))
  
  # #0.5*log(|Q_ac|)
  #reo <- attr(R_ac, 'pivot')
  
  L_ac <- Matrix::Cholesky(Qstar_ac,perm=TRUE,super=T,LDL = F)
  # mu_hat^T Q_ac_ac mu_hat
  #mu_hat = mustar_ac - Matrix::solve(Qstar_ac,Qstar_aac_v)
  #w <- Matrix::t(mu_hat)%*%Qstar_ac%*%mu_hat #Matrix::solve(Matrix::t(R_ac),mu_hat[reo])
  #lik <- lik - 0.5 * (Matrix::t(w)%*%w)
  lik <- lik - 0.5 *( Matrix::t(mustar_ac)%*%Qstar_ac%*%mustar_ac)
  lik <- lik +  Matrix::t(mustar_ac)%*%Qstar_aac_v
  #w <- Matrix::solve(t(R_ac),Qstar_aac_v)
  w <- Matrix::solve(L_ac,Matrix::solve(L_ac,Qstar_aac_v,system='P'),system='L')
  lik <- lik - 0.5  *Matrix::t(w)%*%w
  
  #lik <- lik + sum(log(diag(R_ac)))
  lik <- lik + Matrix::determinant(L_ac)$mod
  
  
  # y*^T * y
  lik <- lik - sigma_Y^(-2) * 0.5 * input_data$ysTys
  # 1/sigma_Y^m
  lik <- lik - length(input_data$y) * log(sigma_Y)
  
  return(as.vector(lik))
}


#'
#' Computes the loglikelihood (y,b) using the standard likelihood formulation
#' 
#' @param param      - (k x 1) parameter for the model
#' @param input_data - (list) $Bstar
#'                            $ystar
#'                            $T
#'                            $bstar
#'                            $alpha
#'@export                          
likelihood_y_gb2<-function(param, input_data){
  if(length(param)==7){
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- param[5:6] 
    
    sigma_Y <- exp(param[7])  
  } else if(length(param)==5) {
    kappa <- exp(param[1:2])
    tau   <- exp(param[3:4])
    mu    <- c(0,0)
    sigma_Y <- exp(param[5])
  } else {
    kappa <- exp(param[1])*c(1,1)
    tau   <- exp(param[2])*c(1,1)
    mu    <- c(0,0)
    sigma_Y <- exp(param[3])
  }
  
  
  if(input_data$alpha==2){
    Q_1 <- kappa[1]^4*input_data$C0 + 2*kappa[1]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_1 <-  tau[1]*Q_1
    Q_2 <- kappa[2]^4*input_data$C0 + 2*kappa[2]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
    Q_2 <- tau[2]*Q_2  
  } else if(input_data$alpha==3){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==4){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  } else if(input_data$alpha==5){
    L = kappa[1]^2*input_data$C0 + input_data$G
    Q_1 = tau[1]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
    L = kappa[2]^2*input_data$C0 + input_data$G
    Q_2 = tau[2]*L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L%*%input_data$C0i%*%L
  }
  Sigma_1 <- solve(Q_1)
  Sigma_2 <- solve(Q_2)
  Sigma <- bdiag(Sigma_1,Sigma_2)
  # likelihood A_conX,A_obsX
  #
  Sigma_Y <- input_data$A_obs%*%Sigma%*%Matrix::t(input_data$A_obs)
  Matrix::diag(Sigma_Y) <- Matrix::diag(Sigma_Y) + sigma_Y^2
  Sigma_b <- input_data$A_con%*%Sigma%*%Matrix::t(input_data$A_con)
  Cov_Yb  <- input_data$A_obs%*%Sigma%*%Matrix::t(input_data$A_con)
  
  Sigma_Yb <- rbind( cbind(Sigma_Y, Cov_Yb),
                     cbind(Matrix::t(Cov_Yb), Sigma_b))
  mu_v <- c(rep(mu[1], dim(Q_1)[1]),
            rep(mu[2], dim(Q_2)[1]))
  
  
  mu_Yb <- as.vector(rbind(input_data$A_obs%*%mu_v,input_data$A_con%*%mu_v))
  R <- chol(forceSymmetric(Sigma_Yb))
  w <- Matrix::solve(Matrix::t(R), c(input_data$y, input_data$b)-mu_Yb)
  lik <- -(0.5 * t(w)%*%w) - sum(log(diag(R)))
  
  # - log(pi(b))
  R <- chol(forceSymmetric(Sigma_b))
  w <- Matrix::solve(Matrix::t(R),  input_data$b-input_data$A_con%*%mu_v)
  lik <-  lik +  (0.5 * Matrix::t(w)%*%w) + sum(log(diag(R)))
  return(as.vector(lik))
}


mean_xy2 <- function(param, input_data){
  
  kappa <- exp(param[1:2])
  tau   <- exp(param[3:4])
  mu    <- param[5:6] 
  
  sigma_Y <- exp(param[7])
  Q_1 <- kappa[1]^4*input_data$C0 + 2*kappa[1]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
  Q_1 <-  tau[1]*Q_1
  Q_2 <- kappa[2]^4*input_data$C0 + 2*kappa[2]^2*input_data$G + input_data$G%*%input_data$C0i%*%input_data$G
  Q_2 <- tau[2]*Q_2
  Sigma_1 <- solve(Q_1)
  Sigma_2 <- solve(Q_2)
  Sigma <- bdiag(Sigma_1,Sigma_2)
  # likelihood A_conX,A_obsX
  #
  Sigma_Y <- input_data$A_obs%*%Sigma%*%Matrix::t(input_data$A_obs)
  Matrix::diag(Sigma_Y) <- Matrix::diag(Sigma_Y) + sigma_Y^2
  Sigma_b <- input_data$A_con%*%Sigma%*%Matrix::t(input_data$A_con)
  Cov_Yb  <- input_data$A_obs%*%Sigma%*%Matrix::t(input_data$A_con)
  
  Sigma_Yb <- rbind( cbind(Sigma_Y, Cov_Yb),
                     cbind(Matrix::t(Cov_Yb), Sigma_b))
  mu_v <- c(rep(mu[1], dim(Q_1)[1]),
            rep(mu[2], dim(Q_2)[1]))
  
  
  mu_Yb <- as.vector(rbind(input_data$A_obs%*%mu_v,input_data$A_con%*%mu_v))
  R <- chol(Sigma_Yb)
  Sigma_X_yb <- Sigma%*%cbind(Matrix::t(input_data$A_obs),Matrix::t(input_data$A_con))
  Ex <- Matrix::solve(R,Matrix::solve(Matrix::t(R), c(input_data$y, input_data$b)-mu_Yb))
  Ex <- Sigma_X_yb%*%Ex + mu_v
  return(Ex)
}



cov_deriv <- function(x,y,nu,kappa,sigma){
  if(nu==Inf){
    M1 = (1-kappa^2*y^2)*kappa^2
    M2 = x*y*kappa^4
    M4 = (1-x^2*kappa^2)*kappa^2
    M5 = exp(-0.5*(x^2+y^2)*kappa^2)
    C <- sigma^2*rbind(cbind(M1,M2),cbind(M2,M4))*rbind(cbind(M5,M5),cbind(M5,M5))
  } else {
    h = sqrt(x^2+y^2)
    
    Cv1k = rSPDE::matern.covariance(h,nu=nu-1,kappa=kappa,sigma=sigma)
    #Cv2k = matern.covariance(h,nu=nu-2,kappa=kappa,sigma=sigma)
    Cv2k = (sigma^2 / (2^(nu - 3))) * ((kappa*abs(h))^(nu-2)) * besselK(kappa*abs(h), nu-2)
    Cv2k[h == 0] = sigma^2
    
    M1 = -(kappa^2*gamma(nu-1)/(2*gamma(nu)))*Cv1k + (kappa^4*y^2/(4*gamma(nu)))*Cv2k
    M2 = -(kappa^4*x*y/(4*gamma(nu)))*Cv2k
    M4 = -(kappa^2*gamma(nu-1)/(2*gamma(nu)))*Cv1k + (kappa^4*x^2/(4*gamma(nu)))*Cv2k
    C <- -rbind(cbind(M1,M2),cbind(M2,M4))
  }
  return(C)
}


loglike_cov <- function(sg,input_data){

  sigma=exp(sg[1])
  kappa=exp(sg[2])
  sigma_n=exp(sg[3])
  # regularization
  prior <- - 20* kappa - 2/sigma - 2*log(sigma)
  K <- cov_deriv(input_data$rx,input_data$ry,input_data$nu,kappa,sigma) + sigma_n^2*matlab::eye(length(input_data$y))
  R = chol(K)
  alpha = solve(K,input_data$y)  
  return(0.5*t(input_data$y)%*%alpha + sum(log(diag(R))) + log(2*pi) - prior)
}


mean_cov <- function(param,input_data){
  
  sigma   = exp(param[1])
  kappa   = exp(param[2])
  sigma_n = exp(param[3])
  
  K_cust <- cov_deriv(input_data$rx,input_data$ry,input_data$nu,kappa,sigma) + sigma_n^2*matlab::eye(2*dim(input_data$rx)[1])
  k_cust <- cov_deriv(input_data$rxp,input_data$ryp,input_data$nu,kappa,sigma)
  
  L_cust = chol(K_cust)
  alpha_cust = solve(L_cust,solve(t(L_cust),input_data$y))
  meanCust = t(k_cust)%*%alpha_cust
}


  