matern.bc <- function(x,y,kappa,nu,sigma,bc='none',K = 10,L=1){
  C = 0
  if(bc == 'Dirichlet'){
    for(k in -K:K){
      C = C + matern.covariance(abs(x+2*L*k-y),kappa=kappa,nu=nu,sigma=sigma) -
        matern.covariance(abs(x+2*L*k+y),kappa=kappa,nu=nu,sigma=sigma)
    }
  } else if(bc == 'Neumann'){
    for(k in -K:K){
      C = C + matern.covariance(abs(x+2*L*k-y),kappa=kappa,nu=nu,sigma=sigma) +
        matern.covariance(abs(x+2*L*k+y),kappa=kappa,nu=nu,sigma=sigma)
    }
  } else if(bc == 'Periodic'){
    for(k in -K:K){
      C = C + matern.covariance(abs(x+L*k-y),kappa=kappa,nu=nu,sigma=sigma)
    }
  } else if(bc == 'none'){
    C = matern.covariance(abs(x-y),kappa=kappa,nu=nu,sigma=sigma)
  }
  return(C)
}


matern.derivative <- function(h, kappa, nu, sigma,deriv=1)
{
  if(deriv==1){
    C = h*matern.covariance(h,kappa=kappa,nu=nu-1,sigma=sigma)
    C[h==0] = 0
  } else if (deriv == 2){
    C = matern.covariance(h,kappa=kappa,nu=nu-1,sigma=sigma)+
      h*matern.derivative(h,kappa=kappa,nu=nu-1,sigma=sigma,deriv=1)

  } else {
    C = (deriv-1)*matern.derivative(h,kappa=kappa,nu=nu-1,sigma=sigma,deriv=deriv-2) +
      h*matern.derivative(h,kappa=kappa,nu=nu-1,sigma=sigma,deriv=deriv-1)
  }
  return(-(kappa^2/(2*(nu-1)))*as.matrix(C))
}

fem1d <- function(x,bc)
{
  n = length(x)
  d <- c(Inf, diff(x))
  dm1 = c(d[2:n], Inf)
  G = -bandSparse(n = n, m = n, k = c(-1, 0, 1),
                  diagonals = cbind(1/dm1, -(1/dm1 + 1/d), 1/dm1))
  C = bandSparse(n = n, m = n, k = c(-1, 0, 1),
                 diagonals = cbind(dm1/6, (dm1 + d)/3, d/6))
  if(bc=="Neumann"){
    C[1, 1:2] <- c(d[2], d[2]/2)/3
    C[n, (n - 1):n] <- c(d[n]/2, d[n])/3
  } else if(bc == "Dirichlet"){
    C <- C[2:(n-1),2:(n-1)]
    G <- G[2:(n-1),2:(n-1)]
  } else if(bc == "Periodic"){
    C[1, 1:2] <- c(d[2]+d[n], d[2]/2)/3
    C[1, n-1] <- d[n]/6
    C[n-1, 1] <- d[2]/6
    C = C[1:(n-1),1:(n-1)]
    G[1,1] = 1/dm1[1] + 1/d[n]
    G[1,n-1] = -1/d[n]
    G[n-1,1] = -1/dm1[1]
    G = G[1:(n-1),1:(n-1)]
  }
  return(list(G = G, C = C))
}

matern.neumann <- function(t,s,kappa,nu,T){
  r0 <- function(x){(1+kappa*x)*exp(-kappa*x)}
  M <- matrix(0,nrow=length(t),ncol=length(s))
  C <- exp(-kappa*T)/(2*sinh(kappa*T))
  C1 <- 2*kappa*T/sinh(kappa*T)^2
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      if(nu == 1/2){
        #M[i,j] <- exp(-kappa*abs(t[i]-s[j])) + c(exp(-kappa*t[i]),exp(-kappa*(T-t[i])))%*%matrix(c(1,exp(-kappa*T),exp(-kappa*T),1),nrow=2)%*%c(exp(-kappa*s[j]),exp(-kappa*(T-s[j])))/(1-exp(-2*kappa*T))
        M[i,j] <- (cosh(kappa*(T-abs(t[i]-s[j]))) + cosh(kappa*(t[i]+s[j]-T)))/sinh(kappa*T)
      } else if (nu == 3/2){
        h = t[i]-s[j]
        v = t[i]+s[j]
        M[i,j] <- r0(abs(h)) + C*(r0(h)+r0(-h) +exp(2*kappa*T)*r0(v)+r0(-v)) + C1*cosh(kappa*s[j])*cosh(kappa*t[i])


      }
    }
  }
  return(c(M))
}

matern.neumann.d <- function(t,s,kappa,nu,T){
  r1 <- function(x){-kappa^2*x*exp(-kappa*x)}
  M <- matrix(0,nrow=length(t),ncol=length(s))
  C <- exp(-kappa*T)/(2*sinh(kappa*T))
  C1 <- 2*kappa*T/sinh(kappa*T)^2
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      if(nu == 1/2){

      } else if (nu == 3/2){
        h = -t[i]+s[j]
        v = t[i]+s[j]
        M[i,j] <- sign(h)*r1(abs(h)) +C*(r1(h)-r1(-h)+exp(2*kappa*T)*r1(v)-r1(-v)) + kappa*C1*sinh(kappa*s[j])*cosh(kappa*t[i])
      }
    }
  }
  return(c(M))
}

matern.neumann.d2 <- function(s,t,kappa,nu,T){
  r2 <- function(x){kappa^3*x*exp(-kappa*x)-kappa^2*exp(-kappa*x)}
  M <- matrix(0,nrow=length(t),ncol=length(s))
  C <- exp(-kappa*T)/(2*sinh(kappa*T))
  C1 <- 2*kappa*T/sinh(kappa*T)^2
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      if(nu == 1/2){

      } else if (nu == 3/2){
        h = -t[i]+s[j]
        v = t[i]+s[j]
        M[i,j] <-  -r2(abs(h)) + C*(-r2(h)-r2(-h)+exp(2*kappa*T)*r2(v)+r2(-v)) + C1*kappa^2*sinh(kappa*s[j])*sinh(kappa*t[i])
      }
    }
  }
  return(c(M))
}

matern.n <- function(s,t,kappa,bc = 1){
  if(bc==0){
    C = (1+kappa*abs(s-t))*exp(-kappa*abs(s-t)) + (1+kappa*(s+t))*exp(-kappa*(s+t))
  } else if(bc==1){
    C = (1+kappa*abs(s-t))*exp(-kappa*abs(s-t)) + (1+kappa*(2-s-t))*exp(-kappa*(2-s-t))
  }
  return(C)
}

#d/ds
matern.dn <- function(s,t,kappa,bc = 1){
  if(bc==0){
    C = -kappa^2*(t-s)*exp(-kappa*abs(s-t)) - kappa^2*(s+t)*exp(-kappa*(s+t))
  } else if(bc==1){
    C = -kappa^2*(t-s)*exp(-kappa*abs(s-t)) + kappa^2*(2-s-t)*exp(-kappa*(2-s-t))
  }
  return(C)
}

#d2/dsdt
matern.d2n <- function(s,t,kappa,bc = 1){
  if(bc==0){
    C = kappa^2*(-1+kappa*abs(t-s))*exp(-kappa*abs(s-t)) + kappa^2*(-1+kappa*(s+t))*exp(-kappa*(s+t))
  } else if(bc==1){
    C = kappa^2*(-1+kappa*abs(t-s))*exp(-kappa*abs(s-t)) + kappa^2*(-1+kappa*(2-s-t))*exp(-kappa*(2-s-t))
  }
  return(C)
}
