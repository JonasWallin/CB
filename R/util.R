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
    Ar = rankMatrix(A,method = "qr")
    if(verbose){
      cat(i,": sample", j,", Arank = ", Arank, ", Ar = ",Ar, "OK = ", OK, "\n")
      #cat(OK,"\n")y
    }
    
    
    if(Ar<=Arank){
      OK = setdiff(OK,j)
      A[k,tri] = 0
    } else {
      Arank = Ar
      obs.loc[k,] = u[1]*mesh$loc[tri[1],1:2] + u[2]*mesh$loc[tri[2],1:2] + u[3]*mesh$loc[tri[3],1:2]
      k = k+1
    }
  }
  return(list(A=A,loc=obs.loc))
}
  