# Multi-view CCA
# 
# Author: Daniel V. Samarov
###############################################################################


structgcca <- function(X){
	
	## Number of "views"
	M <- length(X)
	Xcomb <- do.call('cbind', X)
	
	## Get dimensions of each X
	dims <- do.call('rbind', lapply(X, dim))
	nc <- dims[,2]
	
	Cmat <- t(Xcomb) %*% Xcomb
	Dlist <- lapply(X, function(Xi){
				t(Xi) %*% Xi
			})
	
	Dmat <- matrix(0, nrow(Cmat), ncol(Cmat))
	
	k <- 1
	for(i in 1:M){
		ir <- k:(nc[i]+(k-1))
		Dmat[ir,ir] <- Dlist[[i]]
		k <- k+nc[i]
	}
	
	## Left hand inverse
	O <- 1/(M-1)*(Cmat - Dmat)
	Osvd <- matpow(O, 0.5)
	Op <- Osvd$mat
	Oi <- matpow(O, -0.5, Osvd$SVD)$mat
	
	## Get inverse of D matrix
	Dsvd <- matpow(Dmat, -1)
	Di <- Dsvd$mat
	
	## Matrix for running actual decomp
	Smat <- Op %*% Di %*% Op
	
	## Get svd of Smat
	Ssvd <- svd(Smat)
	x
	## Rescale canonical vectors
	Cvec <- Oi %*% Ssvd$v
	
	Pvec1 <- X[[1]] %*% Ssvd$v[1:nc[1], 1:nc[1]]
}

matpow <- function(x, pow, SVD = NULL){
	
	if(is.null(SVD)){
		sx <- svd(x)
	}
	else
		sx <- SVD
	
	U <- sx$u
	V <- sx$v
	D <- diag(sx$d^pow)
	
	out <- U %*% D %*% t(V)
	
	return(list(mat = out, SVD = sx))
}
