Nlm <- function(k,m)
{
	return(sqrt((2*k+1)*(k^2-m^2)/(2*k-1)))
}

SHgradient <- function(x,y,z,Y,Y1,k,m)
{
	G = cbind(-m*y*Y[,k+1-m] - Nlm(k,m)*x*z*Y1[,k+1+m] + k*x*z^2*Y[,k+1+m],
			   m*x*Y[,k+1-m] - Nlm(k,m)*y*z*Y1[,k+1+m] + k*y*z^2*Y[,k+1+m],
			   Nlm(k,m)*Y1[,k+1+m] - k*z*Y[,k+1+m])
	G[,1] = G[,1]/(1-z^2)
	G[,2] = G[,2]/(1-z^2)
	return(t(G))
}

vcrossp <- function(a, b) 
{ 
   result <- matrix(NA, 3,ncol(a)) 
   result[1,] <- a[2,]*b[3,] - a[3,]*b[2,] 
   result[2,] <- a[3,]*b[1,] - a[1,]*b[3,] 
   result[3,] <- a[1,]*b[2,] - a[2,]*b[1,] 
   return(result) 
} 

SHcross <- function(x,y,z,Y,Y1,k,m)
{
	H = cbind(k*z*y*Y[,k+1+m]- Nlm(k,m)*y*Y1[,k+1+m]+m*x*z*Y[,k+1-m],
			  -k*z*x*Y[,k+1+m]+ Nlm(k,m)*x*Y1[,k+1+m]+m*y*z*Y[,k+1-m],
			   -m*Y[,k+1-m])
	H[,1] = H[,1]/(1-z^2)
	H[,2] = H[,2]/(1-z^2)
	return(t(H))
}

vector.basis <- function(mesh,ord,fix.ord=FALSE,FVc=FALSE,rot.inv=FALSE)
{
	P = t(mesh$loc)
	if(FVc){
		FV = mesh$graph$tv
   	 	nF = dim(FV)[1]
   	 	Ptri = sapply(1:nF, function(f) rowMeans(P[,FV[f,]]))
    	P = sapply(1:nF, function(f) Ptri[,f]/sqrt(t(Ptri[,f])%*%Ptri[,f]))
	}

	x = P[1,]
	y = P[2,]
	z = P[3,]
	i1= i2=1
	n = dim(P)[2]

	#if(rot.inv){
	#	Sph = inla.fmesher.smorg(t(as.array(P)), mesh$graph$tv, sph0 = ord)$sph
	#} else {
	Sph = inla.fmesher.smorg(t(as.array(P)), mesh$graph$tv, sph = ord)$sph
    #}    
	B=list()
	for(i in 1:3){
		if(rot.inv){
			B[[i]] =  Matrix(0,n,2*ord)
		} else {
			B[[i]] =  Matrix(0,n,2*dim(Sph)[2]-1)
		}
	}

	if(fix.ord){
		ord1 = ord-1
	} else {
		ord1 = ord
	}
	for(kk in 1:ord1) {
		if(0){#rot.inv){
			Y = Sph[,kk]
			Y1 = Sph[,max(kk-1,1)]
			ind.set = 0
		} else {
			ind.set = 0#-kk:kk
		    Y = Sph[,(sum(seq(1,1+2*(kk-1),by=2))+1):sum(seq(1,1+2*kk,by=2))]
   		 	if(kk>1){
   		 Y1 = Sph[,(sum(seq(1,1+2*(kk-2),by=2))+1):sum(seq(1,1+2*(kk-1),by=2))]
			} else {
				Y1 = Sph[,1]
			}
		}
		Y1 = cBind(rep(0,n), Y1, rep(0,n))		

	    for(mm in ind.set) {
   	 		G = SHgradient(x,y,z,Y,Y1,kk,mm)
   	     	for(j in 1:3){
            	B[[j]][,i2]=G[j,]
        	}
     	   if(kk<ord-1 && fix.ord) {
        	    H = vcrossp(G,P)
           		for(j in 1:3) {  
                	B[[j]][,i2+1] = H[j,]
            	}
            	i2=i2+1
        	} else if (!fix.ord){
        		H = SHcross(x,y,z,Y,Y1,kk,mm)
            	for(j in 1:3) {  
               		B[[j]][,i2+1] = H[j,]
            	}
            	i2=i2+1
        	}
        	i2=i2+1
        	i1=i1+1
    	}
	}

	for(j in 1:3){
    	B[[j]][is.na(B[[j]])==1]=0
    	B[[j]] = t(B[[j]])
	}
	return(list(B=B,P=P))
}

create.vector.field <- function(B,theta.vec)
{
	return(cbind(colSums(diag(theta.vec)%*%B[[1]]),
				 colSums(diag(theta.vec)%*%B[[2]]),
				 colSums(diag(theta.vec)%*%B[[3]])))
}

vector.basis.2d <- function(mesh,FVc=FALSE)
{
	if(FVc){
		n = dim(mesh$graph$tv)[1]
		FV = mesh$graph$tv
    	nF = dim(FV)[1]
    	Ptri = sapply(1:nF, function(f) rowMeans(P[,FV[f,]]))
    	P = sapply(1:nF, function(f) Ptri[,f]/sqrt(t(Ptri[,f])%*%Ptri[,f]))
	} else {
		n = dim(mesh$loc)[1]
		P = mesh$loc
	}
	B=list()
	for(i in 1:3) {
		B[[i]] =  Matrix(0,3,n)
		B[[i]][i,] = 1
	}
	return(list(B=B,P=P))
}

plot.nested.cov <- function(ind,Q,H,proj, main="")
{
	
	Sigma = (H %*% solve(Q,Matrix(H[ind,],dim(Q)[1],1)))[,1]
	image.plot(proj$x,proj$y,inla.mesh.project(proj,field=Sigma),
				xlab="",ylab="", main = main)
}

plot.vector.field <- function(P,vec,vec.scale, vec.length,S2=FALSE, thin, main="",col="black",add=FALSE)
{
	if(dim(P)[1]< dim(P)[2]){
		P = t(P)
	}
	if(!missing(thin)){
	    ind = seq(1,dim(P)[1],by=thin) 
    	P=P[ind,]
    	vec=vec[ind,]
	}
	if(S2){
    P2 = P + vec.scale*vec
    P2 = P2/sqrt(rowSums(P2^2))
    th1 = asin(P[,3])
    ph1 = atan2(P[,2],P[,1])
    th2 = asin(P2[,3])
    ph2 = atan2(P2[,2],P2[,1])  
    p1 <- mapproject(th1, ph1, projection="mollweide")
    p2 <- mapproject(th2, ph2, projection="mollweide")

    l1 = cbind(p1$x,p1$y)
    l2 = cbind(p2$x,p2$y)
    ld = l2 - l1
    ind1 = (rowSums(ld^2) < 0.004)
    ind2 = (rowSums(ld^2) > 1e-8)
    ind = ind1&ind2
    plot(NULL, type = "n", xlim=c(min(p1$x),max(p1$x)), 
        ylim=c(min(p1$y),max(p1$y)))
    arrows(p1$x[ind], p1$y[ind], p2$x[ind], p2$y[ind],length=vec.length)
	} else {
	  P2 = P + vec.scale*vec 
    ind = colSums((P2-P)^2) > 1e-8
    if(!add){
      plot(NULL, type = "n", xlim=c(min(P[,1]),max(P[,1])), 			
           ylim=c(min(P[,2]),max(P[,2])), main = main)  
    }
		
		arrows(P[ind,1], P[ind,2], P2[ind,1], P2[ind,2],length=vec.length,col=col)
	}
}

make.nested.spde <- function(mesh, spde,B.vec, prior.vec.mu, prior.vec.Q)
{
	n.vec = dim(B.vec[[1]])[1]
	if(missing(prior.vec.mu)) prior.vec.mu = rep(0,n.vec)
	if(missing(prior.vec.Q)){ 
		prior.vec.Q = sparseMatrix(i=1:n.vec,j=1:n.vec,x = rep(0.01,n.vec),
									dims = c(n.vec,n.vec)) 
	}
	fem = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2, gradients=TRUE)
	
	n = dim(B.vec[[1]])[2]
	Bvec = list()
	for(i in 1:n.vec){
	B.vec[[i]]=sparseMatrix(i=1:n,j=1:n,x=B.vec[[1]][i,],dims=c(n,n))%*%fem$dx + 
			sparseMatrix(i=1:n,j=1:n,x=B.vec[[2]][i,],dims=c(n,n))%*%fem$dy +
			sparseMatrix(i=1:n,j=1:n,x=B.vec[[3]][i,],dims=c(n,n))%*%fem$dz 
	}
	
	prior.vec = list(mu = prior.vec.mu, Q = prior.vec.Q)
	nested.spde = list(spde = spde, 
					   model = "nested spde",
					   B.vec = B.vec,
					   n.vec = n.vec,
					   prior.vec = prior.vec)
	return(nested.spde)
}

nested.spde.H <- function(nested.spde,theta.vec)
{
	n = dim(nested.spde$spde$param.inla$M0)[1]
	Hl = 0
	if(!is.na(theta.vec[1])){
		for(i in 1:length(theta.vec)){
			Hl = Hl + theta.vec[i]* nested.spde$B.vec[[i]]
		}
	}
	H = Hl + sparseMatrix(i=1:n,j=1:n,x=rep(1,n),dims=c(n,n))
	return(H)
}

nested.spde.sample <- function(nested.spde,theta.spde,theta.vec){
	Q = inla.spde2.precision(nested.spde$spde, theta=theta.spde)
	x = inla.qsample(Q=Q)
	y = nested.spde.H(nested.spde,theta.vec)%*%x
	return(list(x=x,y=y))
}

nested.spde.stack <- function(Y,A,B,sigma2.hp, mean.mu, mean.Q)
{
	if(missing(sigma2.hp)) sigma2.hp = c(2, 500)
	n.cov = dim(B)[2]
	if(is.null(n.cov)) n.cov = 1
	if(missing(mean.mu)) mean.mu = rep(0,n.cov)
	if(missing(mean.Q)){ 
		mean.Q = sparseMatrix(i=1:n.cov, j=1:n.cov, 
							x= rep(0.01,n.cov), dims = c(n.cov,n.cov))
	}
 	
	return(list(Y=Y, A=A, B=B, n.obs = length(Y), n.mu = n.cov,
				 sigma2 = list(a = sigma2.hp[1], b= sigma2.hp[2]), 
				 mu=list(mu = mean.mu, Q = mean.Q), 
				 AtA = t(A) %*% A, YtY = t(Y) %*% Y, 
				 AtB = t(A) %*% B, AtY = t(A) %*% Y, 
				 BtY = t(B) %*% Y, BtB = t(B) %*% B))
}
