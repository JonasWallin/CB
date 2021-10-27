#' @title Sampling of locations in triangulations
#'
#' @description Create observation locations in triangulation by sampling
#' uniformly over triangles and then uniformly within each triangle.
#' The locations are sampled so that the rank of the observation matrix is n.obs
#'
#' @param mesh the mesh.
#' @param n.obs the total number of observations needed.
#' @return a list with the observation matrix (A),
#'     the observation locations (loc), and the number observations per triangle (n.obs)
#' @export
#'
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

#' @title Sampling of locations in triangulations
#'
#' @description Create observation locations in triangulation by sampling
#' uniformly over triangles and then uniformly within each triangle.
#' Only one observation per triangle is allowed.
#'
#' @param mesh the mesh.
#' @param n.obs the total number of observations needed.
#' @return a list with the observation matrix (A),
#'     the observation locations (loc), and the number observations per triangle (n.obs)
#' @export
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


#' @title Mesh construction for example on GP regression
#' @description Function used to construct a mesh for the constrained GP regression
#'
#' @param grid_ (2 x 1) start end location of grid
#' @param expan (double) expand the boundary to remove boundary effect
#' @param n     (int) number of mesh points
#' @param const.by - (int) often to evalute gradient
#' @param soft_dirchlet (bool) use to remove boundary effects
#' @param soft_derivative (bool) use soft constraint
#' @return  input_data (list) for likelhiood or mean construction
#' @export
build_mesh_divergence <- function(grid_,
                                  expand,
                                  n,
                                  const.by=3,
                                  soft_dirchlet=TRUE,
                                  soft_derivatieve=TRUE,
                                  alpha)
{
    if(2*expand > diff(grid_)){
        #large extension, use graded mesh so that we have at most n/2 nodes in the extension
        x.int = seq(from=grid_[1],to=grid_[2],length.out=ceil(n/2))
        h = x.int[2]-x.int[1]
        a = (expand/h-n/4)/(n/4+1)
        tmp1 <- (1+seq(from=0,to=a,length.out = floor(n/4)))
        tmp <- c(fliplr(tmp1), rep(1,n-1-2*length(tmp1)), tmp1)
        h <- (diff(grid_) + 2*expand)/sum(tmp)
        tmp <- h*c(0,tmp)
        x_loc <- grid_[1]-expand + cumsum(tmp)
    } else {
        x_loc = seq(from = grid_[1]-expand, to = grid_[2]+expand,length.out = n)
    }

    mesh <- INLA::inla.mesh.1d(x_loc)
    P = t(mesh$loc)
    Sph =INLA::inla.mesh.fem(mesh, order = 1)
    input_data <- list()

    h = diff(x_loc)
    C1 = Sph$c1
    G1 = Sph$g1

    i.l = 2:mesh$n
    i.r = 1:(mesh$n - 1)
    i.0 = 1:mesh$n
    j.l = 1:(mesh$n - 1)
    j.r = 2:mesh$n
    j.0 = 1:mesh$n
    v.l <- -0.5*h^2
    v.r <- 0.5*h^2
    v.0 <- rep(0,mesh$n)
    v.0[1] = -0.5*h[1]^2
    v.0[mesh$n] = -0.5*h[mesh$n-1]^2
    H1 = Matrix::sparseMatrix(i = c(i.l, i.r, i.0), j = c(j.l, j.r, j.0), x = c(v.l, v.r, v.0), dims = c(mesh$n, mesh$n))
    C1 <- C1[2:(n-1),2:(n-1)]
    G1 <- G1[2:(n-1),2:(n-1)]
    H1 <- H1[2:(n-1),2:(n-1)]

    n <- n-2

    C <- kronecker(C1,C1)

    G <- kronecker(G1,C1) + kronecker(C1,G1)

    C0 = diag(rowSums(C1))
    C0i = diag(1/rowSums(C1))
    Hx <- kronecker(H1,C1)
    Hy <- kronecker(C1,H1)


    C0 = Matrix::sparseMatrix(i=1:n^2,j=1:n^2,x=rowSums(C),dims=c(n^2,n^2))
    C0i = Matrix::sparseMatrix(i=1:n^2,j=1:n^2,x=1/rowSums(C),dims=c(n^2,n^2))


    I = Matrix::Matrix(0,n,n, sparse=TRUE)
    index_grid <- seq(2,n-1,by=const.by)
    in_grid <- 1:n#which(grid_[1]-expand/2<x_loc[2:(n+1)] & x_loc[2:(n+1)] < grid_[2]+expand/2)
    index_grid <- intersect(index_grid,in_grid)
    I[index_grid, index_grid] <- 1
    index <- which(I>0)
    A_con <- cbind(Hx[index,],Hy[index, ])
    index_0 <- which(I==0)
    A_con_fake <- cbind(Hx[index_0,],Hy[index_0, ])

    SVD_Build <- c_basis2(A_con)
    b <- rep(0, dim(A_con)[1])
    Tmat <- SVD_Build$T
    a <- 1:(dim(A_con)[1])
    ac <- (dim(A_con)[1]+1):(2*n^2)
    bstar = rep(0,dim(A_con)[1])

    input_data$x_loc <- x_loc
    input_data$n <- n
    input_data$mesh <- mesh
    input_data$alpha <- alpha
    input_data$C0 <- C0
    input_data$G  <- G
    input_data$C0i <- C0i
    input_data$Tmat <- Tmat
    input_data$ac <- ac
    input_data$a <- a
    input_data$A_con <- A_con
    input_data$b     <- rep(0, dim(A_con)[1])
    input_data$bstar <-  (t(SVD_Build$U)%*%b)/SVD_Build$S

    # fake obs loc
    if(soft_dirchlet || soft_derivatieve)
    {
        grid_fake <- c(grid_[1]-expand/2,grid_[2]+expand/2)
        ind <- which(grid_fake[1]<x_loc& x_loc< grid_fake[2])
        A1 <- INLA::inla.mesh.1d.A(mesh,x_loc[ind])
        A1 <- A1[,2:(n+1)]
        l1 <- length(x_loc[ind])
        p1 <- x_loc[which(abs(x_loc-grid_fake[1])==min(abs(x_loc-grid_fake[1])))]
        p2 <- x_loc[which(abs(x_loc-grid_fake[2])==min(abs(x_loc-grid_fake[2])))]
        A2 <- INLA::inla.mesh.1d.A(mesh,rep(p1,l1))
        A3 <- INLA::inla.mesh.1d.A(mesh,rep(p2,l1))
        A2 <- A2[,2:(n+1)]
        A3 <- A3[,2:(n+1)]
        A4 <- rbind(A1,A1,A2,A3)
        A5 <- rbind(A2,A3,A1,A1)
        A_fake <- t(KhatriRao(t(A4),t(A5)))
        A_fake <- bdiag(A_fake,A_fake)
        xy_fake <-rbind(cbind(x_loc[ind], p1),
                        cbind(x_loc[ind], p2),
                        cbind(p1, x_loc[ind]),
                        cbind(p2, x_loc[ind]))
        B_fake <- A_fake%*%t(Tmat)
        B_fake_ac <- B_fake[, ac]
        input_data$BtB_fake <- 0*(t(B_fake_ac)%*%B_fake_ac)
        if(soft_dirchlet)
            input_data$BtB_fake   <- 10*(t(B_fake_ac)%*%B_fake_ac)

        B_fake2 <- A_con_fake%*%t(Tmat)
        B_fake2_ac <- B_fake2[,ac]
        if(soft_derivatieve)
            input_data$BtB_fake   <- input_data$BtB_fake + 10^12*(t(B_fake2_ac)%*%B_fake2_ac)


    }

    return(input_data)
}
#' @title Add observations to divergence example
#' @description Utility function to add observation for the example on constrained GP regression
#' @param Y          - data
#' @param obs.loc    - location of observations
#' @param input_data - output from build_mesh_divergence
#' @export
add_observations <-function(Y, obs.loc,input_data ){

    Ax = INLA::inla.mesh.1d.A(input_data$mesh,obs.loc[,1])
    Ay = INLA::inla.mesh.1d.A(input_data$mesh,obs.loc[,2])
    Ax <- Ax[,2:(input_data$n+1)]
    Ay <- Ay[,2:(input_data$n+1)]

    Aobs <- t(KhatriRao(t(Ax),t(Ay)))
    A_obs <- bdiag(Aobs,Aobs)


    input_data$A_obs <- A_obs
    input_data$y     <- Y
    B <- A_obs%*%t(input_data$Tmat)
    Bac <- B[, input_data$ac]
    input_data$BtB   <- t(Bac)%*%Bac
    input_data$BtY   <- t(Bac)%*%(Y - B[,input_data$a]%*%input_data$bstar)
    input_data$ysTys <- t((Y - B[,input_data$a]%*%input_data$bstar))%*%((Y - B[,input_data$a]%*%input_data$bstar))


    return(input_data)
}
