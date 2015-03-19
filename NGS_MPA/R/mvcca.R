# Multi-view CCA
# 
# Author: Daniel V. Samarov
###############################################################################


clusCCA <- function(X, group){
	
	## Number of "views"
	M <- length(X)
	Xcomb <- do.call('cbind', X)
	Cmat <- crossGroupCov(Xcomb, group)
	
	## Get dimensions of each X
	dims <- do.call('rbind', lapply(X, dim))
	nc <- dims[,2]
	
	Dlist <- lapply(X, function(Xi){
				groupCov(Xi, group)
			})
	
	Dmat <- matrix(0, nrow(Cmat), ncol(Cmat))
	
	k <- 1
	for(i in 1:M){
		ir <- k:(nc[i]+(k-1))
		Dmat[ir,ir] <- Dlist[[i]]
		k <- k+nc[i]
	}
	
	## Left hand inverse
	Dsvd <- matpow(Dmat, -0.5)
	Di <- Dsvd$mat
	
	CDmat <- 1/(M-1)*(Cmat - Dmat)
	
	## Matrix for running actual decomp
	Smat <- Di %*% CDmat %*% Di
	
	## Get svd of Smat
	Ssvd <- svd(Smat)
	
	## Rescale canonical vectors
	Cvec <- Di %*% Ssvd$v
	
	return(Cvec)
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

## Cross covariance between clusters
crossGroupCov <- function(x, group){
	
	spX <- split(as.data.frame(x), group)
	ncores <- detectCores()
	cl <- makeCluster(ncores)
	browser()
	CxyList <- lapply(1:length(spX), function(i){
				allPairsCov(as.matrix(spX[[i]]), cl)
			})
	Cxy <- Reduce('+', CxyList)
	stopCluster(cl)
	return(Cxy)
}

allPairsCov <- function(x, cl){
	
	n <- nrow(x)
	nc <- ncol(x)
	browser()
	clusterExport(cl, c('x', 'n', 'nc'))
	tm <- system.time({CxyList <- parLapply(cl, 1:n, fun = function(i){
				xrep <- matrix(rep(x[i,], times = n), n, nc,
						byrow = TRUE)
				cxy <- t(xrep) %*% x
				cxy
			})})
	Cxy <- Reduce('+', CxyList)
	
	return(Cxy)
}

groupCov <- function(x, group){
	
	spX <- split(as.data.frame(x), group)
	covs <- lapply(spX, function(u){
				um <- as.matrix(u)
				nrow(um) * t(um)%*%um
			})
	gcov <- (1/M) * Reduce('+', covs)
	gcov <- nrow(spX[[1]])*covs[[1]]
	
	return(gcov)
}
