# Multi-view CCA
# 
# Author: Daniel V. Samarov
###############################################################################

require(parallel)
require(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp('src/allPairsCov.cpp')
clusCCA <- function(X, group, reg = 0, 
		rnk = NULL, reduce = TRUE, frac = 0.2){
	
	## Number of "views"
	M <- length(X)
	Xcomb <- do.call('cbind', X)
	Cmat <- crossGroupCov(Xcomb, group, frac = frac, 
			reduce = reduce)
	
	## Get dimensions of each X
	dims <- do.call('rbind', lapply(X, dim))
	nc <- dims[,2]
	
	Dlist <- lapply(X, function(Xi){
				groupCov(Xi, group)
			})
	
	Dmat <- matrix(0, ncol(Xcomb), ncol(Xcomb))
	
	k <- 1
	for(i in 1:M){
		ir <- k:(nc[i]+(k-1))
		Dmat[ir,ir] <- Dlist[[i]]
		k <- k+nc[i]
	}
	
	## Left hand inverse
	Dsvd <- matpow(Dmat, -0.5, rnk, reg)
	Di <- Dsvd$mat
	
	CmatPD <- matpow(Cmat, 1, rnk)$mat
	CDmat <- 1/(M-1)*(CmatPD - Dmat)
	
	## Matrix for running actual decomp
	Smat <- Di %*% CDmat %*% Di
	
	## Get svd of Smat
	Ssvd <- svd(Smat)
	
	## Rescale canonical vectors
	Cvec <- Di %*% Ssvd$v
	colnames(Cvec) <- colnames(Xcomb)
	
	return(Cvec)
}

clusCCA.pen <- function(X, group, reg = 0, 
		rnk = NULL){
	
	## Number of "views"
	M <- length(X)
	Xcomb <- do.call('cbind', X)
	Cmat <- t(Xcomb) %*% Xcomb
	
	## Get dimensions of each X
	dims <- do.call('rbind', lapply(X, dim))
	nc <- dims[,2]
	
	Dlist <- lapply(X, function(Xi){
				groupSS(Xi, group)
			})
	
	Dmat <- matrix(0, ncol(Xcomb), ncol(Xcomb))
	
	k <- 1
	for(i in 1:M){
		ir <- k:(nc[i]+(k-1))
		Dmat[ir,ir] <- t(X[[i]])%*%X[[i]]+reg*Dlist[[i]]
		k <- k+nc[i]
	}
	
	## Left hand inverse
	Dsvd <- matpow(Dmat, -0.5, rnk, 0)
	Di <- Dsvd$mat
	
	CmatPD <- matpow(Cmat, 1, rnk)$mat
	CDmat <- 1/(M-1)*(CmatPD - Dmat)
	
	## Matrix for running actual decomp
	Smat <- Di %*% CDmat %*% Di
	
	## Get svd of Smat
	Ssvd <- svd(Smat)
	
	## Rescale canonical vectors
	Cvec <- Di %*% Ssvd$v
	colnames(Cvec) <- colnames(Xcomb)
	
	return(Cvec)
}


matpow <- function(x, pow, rnk = NULL, reg = 0, SVD = NULL){
	
	if(is.null(SVD)){
		sx <- svd(x)
	}
	else
		sx <- SVD
	
	if(is.null(rnk))
		rnk <- ncol(x)
	
	U <- sx$u
	V <- sx$v
	val <- (sx$d+reg)^pow
	val[-c(1:rnk)] <- 0
	D <- diag(val)
	
	if(is.null(rnk))
		rnk <- ncol(x)
	
	out <- U %*% D %*% t(V)
	
	return(list(mat = out, SVD = sx))
}

## Cross covariance between clusters
crossGroupCov <- function(x, group, reduce = TRUE, frac = NULL){
	
	spX <- split(as.data.frame(x), group)
	nc <- ncol(x)
	
	Cxy <- 0
	for(i in 1:length(spX)){
		xi <- as.matrix(spX[[i]])
		if(reduce){
			ns <- ceiling(nrow(xi) * frac)
			xi <- kmeans(xi, ns)$centers
		}
		tmp <- matrix(allPairsCov(xi), nc, nc)
		print(sum((tmp - t(tmp))^2))
		Cxy <- Cxy + tmp
	}
	
	return(Cxy)
}


groupSS <- function(x, group){
	
	spX <- split(as.data.frame(x), group)
	covs <- lapply(spX, function(u){
				um <- scale(as.matrix(u), scale = FALSE)
				t(um)%*%um
			})
	gcov <- Reduce('+', covs)
	
	return(gcov)
}

groupCov <- function(x, group){
	
	spX <- split(as.data.frame(x), group)
	covs <- lapply(spX, function(u){
				um <- scale(as.matrix(u), scale = FALSE)
				nrow(um) * t(um)%*%um
			})
	gcov <- Reduce('+', covs)
	
	Mgrp <- sum(unlist(lapply(spX, nrow))^2)
	return(gcov/Mgrp)
	
}

clusGroups <- function(X, k){
	
	out <- lapply(X, function(x){
				kmeans(x, k)$clus
			})
	return(out)
}

structIndMats <- function(obj){
	
	out <- lapply(obj, function(cl){
				f <- factor(cl)
				model.matrix(~-1+f)
			})
	return(out)
	
}
