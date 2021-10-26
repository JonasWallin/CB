#alg 1
#' build the basis matrix
#' @param A (k x n) the constraint matrix
#' @param T (n x n) th new basis
c_basis1 <- function(A){
  
  n <- dim(A)[2]
  T <- Diagonal(rep(1,n),n=n)
  ind <- abs(Matrix::colSums(abs(A)))>0 #fulhack
  A_id <- A[,ind,drop=FALSE]
  USV <- svd(A_id,nv=dim(A_id)[2])
  T[ind, ind] <- (USV$v)
  T <- cbind(T[,ind,drop=FALSE],T[,ind==FALSE,drop=FALSE]) 
  return(T)
}
#alg 2
#' build the basis matrix
#' @param A (k x n) the constraint matrix
#' @param T (n x n) th new basis
c_basis2 <- function(A){
  n <- dim(A)[2]
  k <- dim(A)[1]
  T <- Diagonal(rep(1,n),n=n)
  index_A <- c(1:k)
  ind <- abs(Matrix::colSums(abs(A)))>0 #fulhack
  A_id <- A[,ind,drop=FALSE]
  TD <- Diagonal(rep(1,sum(ind)),n=sum(ind))
  while(length(index_A)>0){
    index = index_A[1]
    new   = index_A[1]
    while(length(new)>0){
      col_index <- Matrix::colSums(abs(A_id[new,,drop=F]))>0
      full      <- which(Matrix::rowSums(abs(A_id[,col_index,drop=F]))>0)
      new <- setdiff(full, index)
      index <- c(index,new)
    }
    index_col <- Matrix::colSums(abs(A_id[index,,drop=F]))>0
    A_new <- A_id[index, index_col,drop=FALSE]
    USV <- svd(A_new,nv=dim(A_new)[2])
    TD[index_col, index_col] <- USV$v
    index_A <- setdiff(index_A,index)
  }
  T[ind,ind] <- TD
  T <- cbind(T[,ind,drop=FALSE],T[,ind==FALSE,drop=FALSE]) 
  index <- Matrix::colSums(abs((A%*%T)))>10^-10
  T <- cbind(T[,index,drop=FALSE],T[,index==FALSE,drop=FALSE]) 
  return(T)  
}

