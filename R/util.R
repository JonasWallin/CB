#create observation locations in triangulation so that the rank of the observation matrix is n.obs
tri_sample <- function(mesh,n.obs,verbose=FALSE){
  n_tri <- dim(mesh$graph$tv)[1]
  n <- mesh$n
  if(n.obs>n){
    n.obs <- n
  }
  A <- Matrix(data=0,nrow=n.obs,ncol=n,sparse=TRUE)
  obs.loc <- matrix(nrow=n.obs,ncol=2)
  
  OK <- 1:n_tri
  obs.per.tri <- rep(0,n_tri)
  Arank = 0
  k = 1
  i = 0
  while(k <= n.obs){
    i = i+1
    j <- sample(OK,1) #sample one triangle
    tri <- mesh$graph$tv[j,]
    u <- runif(3) # sample uniformly in the triangle
    u = u/sum(u)
    A[k,tri] <- u
    time.start <- Sys.time()
    Ar <- rankMatrix(A,method = "qr")
    if(verbose){
      cat(i,": sample", j,", Arank = ", Arank, ", Ar = ",Ar, "OK = ", OK, "\n")
      #cat(OK,"\n")
    }
    
    
    if(Ar<=Arank){
      OK = setdiff(OK,j)
      A[k,tri] = 0
    } else {
      obs.per.tri[j] = obs.per.tri[j] + 1
      Arank = Ar
      obs.loc[k,] = u[1]*mesh$loc[tri[1],1:2] + u[2]*mesh$loc[tri[2],1:2] + u[3]*mesh$loc[tri[3],1:2]
      k = k+1
      #if(obs.per.tri[j]==1){
        OK = setdiff(OK,j)  
      #}
    }
  }
  return(list(A=A,loc=obs.loc,n.obs = obs.per.tri))
}

# Create observation locations in triangulation by sampling uniformly over triangles and then uniformly within each triangle. 
# Only one observation per triangle is allowed. 
tri_sample_simple <- function(mesh,n.obs,verbose=FALSE){
  n_tri <- dim(mesh$graph$tv)[1]
  n <- mesh$n
  if(n.obs>n){
    n.obs <- n
  }
  A <- Matrix(data=0,nrow=n.obs,ncol=n,sparse=TRUE)
  obs.loc <- matrix(nrow=n.obs,ncol=2)
  
  OK <- 1:n_tri
  obs.per.tri <- rep(0,n_tri)
  k = 1
  while(k <= n.obs){
    j <- sample(OK,1) #sample one triangle
    tri <- mesh$graph$tv[j,]
    u <- runif(3) # sample uniformly in the triangle
    u = u/sum(u)
    A[k,tri] <- u
    obs.per.tri[j] = obs.per.tri[j] + 1
    obs.loc[k,] = u[1]*mesh$loc[tri[1],1:2] + u[2]*mesh$loc[tri[2],1:2] + u[3]*mesh$loc[tri[3],1:2]
    k = k+1
    OK = setdiff(OK,j)  
  }
  return(list(A=A,loc=obs.loc,n.obs = obs.per.tri))
}
