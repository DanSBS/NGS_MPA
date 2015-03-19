# TODO: Add comment
# 
# Author: samarov
###############################################################################


clusCCA <- function (datasets, reg = 0, group) 
{
	mat <- datasets
	m <- length(mat)
	para <- reg
	if (para < 0 || para > 1) 
		stop(" Value of 'reg' should be greater than zero and less than one")
	drop <- function(mat, parameter) {
		mat <- as.matrix(mat)
		p <- parameter
		fac <- svd(mat)
		len <- length(fac$d)
		ord <- order(fac$d)
		ord_eig <- fac$d[ord]
		portion <- cumsum(ord_eig^2)/sum(ord_eig^2)
		ind <- which(portion < p)
		num <- len - length(ind)
		fac$v <- as.matrix(fac$v[, 1:num])
		fac$u <- as.matrix(fac$u[, 1:num])
		mat_res <- (fac$v %*% sqrt(diag(1/fac$d[1:num], nrow = length(fac$d[1:num]))) %*% 
					t(fac$u))
		w <- list(mat_res = mat_res, num = num)
		return(w)
	}
	covm <- array(list(), m)
	fea <- c(rep(0, m))
	mean_m <- array(list(), m)
	
	## This is where the cluster CCA code differs
	for (i in 1:m) {
		mean_m[[i]] <- apply(mat[[i]], 2, mean)
#		mat[[i]] <- apply(mat[[i]], 2, function(x) {
#					x - mean(x)
#				})
		covm[[i]] <- groupCov(mat[[i]],group)
		fea[i] <- ncol(mat[[i]])
	}
	whiten_mat <- array(list(), m)
	cov_whiten <- array(list(), m)
	white <- array(list(), m)
	d <- c(rep(0, m))
	if (para > 0) {
		for (i in 1:m) {
			dummy <- drop(covm[[i]], para)
			whiten_mat[[i]] <- mat[[i]] %*% dummy$mat_res
			d[i] <- dummy$num
			white[[i]] <- dummy$mat_res
		}
	}
	else {
		for (i in 1:m) {
			whiten_mat[[i]] <- mat[[i]] %*% SqrtInvMat(covm[[i]])
			d[i] <- ncol(mat[[i]])
			white[[i]] <- SqrtInvMat(covm[[i]])
		}
	}
	if (sum(d) == sum(fea)) {
		print("Normal gCCA")
	}
	else {
		print("Regularized gCCA")
	}
	a <- concatenate(whiten_mat)
	z <- cov(a)
	eig <- svd(z)
	proj_data <- a %*% eig$v
	eig_wh <- array(list(), m)
	if (para > 0) {
		j <- 0
		for (i in 1:m) {
			dummy <- drop(covm[[i]], para)
			eig_wh[[i]] <- dummy$mat_res %*% eig$v[(j + 1):(j + 
								fea[i]), ]
			j <- j + fea[i]
		}
	}
	else {
		j <- 0
		for (i in 1:m) {
			eig_wh[[i]] <- SqrtInvMat(covm[[i]]) %*% eig$v[(j + 
								1):(j + fea[i]), ]
			j <- j + fea[i]
		}
	}
	gencorr <- list(eigval = eig$d, eigvecs = eig_wh, proj = proj_data, 
			meanvec = mean_m, white = white)
	return(gencorr)
}


groupCov <- function(x, group){
	
	covs <- lapply(split(as.data.frame(x), group), cov)
	gcov <- covs[[1]]
	M <- length(covs)
	for(i in 2:M){
		gcov <- gcov + covs[[i]]
	}
	
	gcov <- (1/M)*gcov
	
	return(gcov)
}