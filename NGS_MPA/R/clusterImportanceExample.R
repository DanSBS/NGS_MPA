# Testing out some concepts
# 
# Daniel V. Samarov
###############################################################################
require(genlasso)
require(ggplot2)
require(MASS)
source("C:\\Users\\samarov\\git\\gsmethods\\gsmethods\\R\\smtlasso_v2.R")
source("C:/eclipse/workspace372/CARS Analysis/Simulation Study/Cell Simulation/R/vca.R")
dataDir <- 'C:/eclipse/workspace/OMI/Dye Mixture Analysis/2 Dye Mix/output/cubes.RData'
cubes <- load(dataDir)
#c1 <- mvrnorm(200, c(-1,0), matrix(c(.1,.25,.25,.9),2,2))
#c2 <- mvrnorm(200, c(1, 0), matrix(c(.1,.25,.25,.9),2,2))
#c3 <- mvrnorm(200, c(-1,-3), matrix(c(.1,.25,.25,.9),2,2))
#x <- rbind(c1,c2,c3)


k <- 5
ngam <- 50
max.iter <- 1000
eps <- 1e-3
lambda <- .1
clus <- FALSE

x0 <- matrix(c(cube2),
		nrow(cube2) * ncol(cube2))
xb <- baseline(x0, method = 'als')
x1 <- t(apply(x0, 1, function(u){
					locpoly(1:61, u, drv = 0,
							degree = 1, bandwidth = 1,
							gridsize = 61)$y
				}))
xb <- t(apply(x1, 1, function(u){
					baseline(rbind(u))@corrected		
				}))
x <- apply(x0, 2, function(u){
			c(image.smooth(matrix(u,128,128), theta = 1)$z)
		})
#x <- t(apply(x1, 1, function(u) (u - min(u))/diff(range(u))))

V <- cbind()
for(i in 2:k){
	V <- cbind(V, vca(x, i, center = FALSE)$Ae)
}

cV <- cor(V)
ind <- c()
cIndx <- rep(TRUE, ncol(V))
for(i in 1:ncol(V)){
	iold <- cIndx[i]
	cIndx[abs(cV[i,]) > .99] <- FALSE
	if(iold)
		ind <- c(ind, i)
}
V <- V[,ind]
#V <- (V - min(V))/diff(range(V)) * max(x)
V <- apply(V, 2, function(u) (u - min(u))/diff(range(u)))
#V <- (V - min(V))/diff(range(V))

V <- vca(x, k, center = TRUE)$Ae
mfcTV <- mfClus(x, k, lambda, alpha = .1, Vstart = V, 
		reduce = TRUE, rsize = 1000, iter = 5, 
		clus = FALSE, normU = FALSE, pos = TRUE)
mfcGrp <- mfClus(x, k, lambda, Vstart = V, 
		reduce = TRUE, rsize = 1000, iter = 5,
		type = 'group', lambdag = .05, lambda = .01)
Vscl <- sinh(t(scale(t(mfcGrp$V), center = FALSE)) * scl[[1]] + cent[[1]])

k <- ncol(mfc$V)
o <- apply(mfc$V, 1, function(v) length(unique(round(v,6)))) 
write.csv(Vscl[order(o, decreasing = TRUE),], 
		file = 'output/ILL250_scaled_v3.csv')


pdf('plots/clusterImportance_ILL250_ROC_v3.pdf',
		height = 8, width = 11)
clusterHeatMap(t(scale(t(mfc$V), center = FALSE)))
dev.off()


mfClus <- function(x, k, 
		lambda, alpha = 0.001, 
		Vstart = NULL, ngam = 100, max.iter = 1000,
		eps = 1e-3, iter = 3,
		reduce = FALSE, rsize = 2000,
		type = 'tv', lambdag = 1, clus = TRUE, normU = FALSE,
		pos = FALSE){
	
	## POS is super hacky right now
	## Get initial starting values as centers 
	## from kmeans
	if(is.null(Vstart))
		V <- t(kmeans(x,k)$cent)
	else{
		V <- Vstart
		k <- ncol(V)
	}
	
	if(reduce){
		km <- kmeans(x, rsize)
		s <- km$clus
		Xred <- km$cen
	}
	else
		Xred <- x
	
	oldDiffs <- 1e10
	
	for(i in 1:iter){
		print(paste('Iteration:', i))
		
		print('=================RUNNING SGL=====================')
		U <- t(smtlasso(V, t(x), NULL, lambda = .01, 
						rho = 0, ngam = ngam,
						max.iter = max.iter,
						parallel = FALSE,
						eps = eps, pos = 1, 
						method='sgl', 
						intercept=FALSE, clus = clus,
						singled=1))
		
		## Remove 0 columns
		chk0 <- colSums(U) != 0
		U <- U[,chk0]
		if(normU){
			rU <- rowSums(U)
			rU[rU == 0] <- Inf
			U <- U/rU
		}
		
		k <- ncol(U)
		print(paste('New # clusters:', k))
		
		if(k < 2)
			stop('Group sparse parameters too large')
		## Check type of penalty structure
		if(type == 'group'){
			V <- t(smtlasso(U, x, NULL, lambda = lambdag, 
							rho = 0, ngam = ngam,
							max.iter = max.iter,
							parallel = FALSE,
							eps = eps, pos = 0, 
							method='sgl', 
							intercept=FALSE,
							singled=1, clus = FALSE))
		}
		if(type == 'tv'){
			D <- diag(k) %x% rep(1, (k - 1))
			oD <- c(rep(1, k) %x% (diag((k - 1))*-1))
			D[D != 1] <- oD
			
			rmInd <- c()
			for(i in 1:(k-1)){
				rmInd <- c(rmInd, (c(i:(k-1))*(k-1) + i))
			}
			D <- -D[-rmInd,]
			
			UD <- rbind(U, 5*D)
			xD <- rbind(x, matrix(0, nrow(D), ncol(x)))
			p <- rbind()
			p2 <- t(solve(t(UD) %*% UD) %*% t(UD) %*% (xD))
			for(i in 1:ncol(xD))
			{
				print(i)
				nn <- nnls(UD, xD[,i])
				p <- rbind(p, coef(nn))
				
			}
			
			if(reduce)
				Ured <- do.call('rbind',lapply(split(as.data.frame(U), s), 
								function(u) colMeans(as.matrix(u))))
			
#				tr <- try({Ured <- kmeans(U, rsize)$cent},TRUE)
			else
				Ured <- U
			
			
			print('=================RUNNING GENLASSO=====================')
			augU <- rbind(Ured, diag(k) * alpha)
			SVD <- svd(augU)
			
			Vlist <- lapply(1:ncol(x), function(i){
						augY <- c(Xred[,i], rep(0, k))
						
						st <- system.time({gl <- genlasso(augY, augU, D,
											maxsteps = 2000,
											approx = TRUE,
											minlam = 1e-5,
											eps = 0.01,
											verbose = FALSE,
											SVD = SVD)})
						gl
					})
			
			print('=================Tuning GENLASSO=====================')
			V <- do.call('rbind', lapply(1:length(Vlist), function(i){
								gl <- Vlist[[i]]
								vi <- gl$beta
								df <- gl$df
								xi <- c(Xred[,i])#, rep(0, k))
								score <- colSums((Ured %*% vi - xi)^2) / df +
										2 * ((nrow(Xred)+0) - df)
								o <- order(score)[1]
								vi[,o]
							}))
#			ll <-  lapply(1:length(Vlist), function(i){
#						gl <- Vlist[[i]]
#						vi <- gl$beta
#						df <- gl$df
#						xi <- c(Xred[,i])#, rep(0, k))
#						score <- colSums((Ured %*% vi - xi)^2) / df +
#								2 * ((nrow(Xred)+0) - df)
#						o <- order(score)[1]
#						vi[,floor(ncol(vi)/2)]
#					})
			if(pos){
				V <- scale(abs(V), center = FALSE)
			}
		}
		
		diffs <- sum((U %*% t(V) - x)^2)/sum(x^2)
		print(paste('Error is:', diffs))
		if(abs(diffs - oldDiffs)/oldDiffs <= 1e-3)
			break
		else
			oldDiffs <- diffs
	}
	
	rownames(V) <- colnames(x)
	
	return(list(U = U, V = V, k = k))
	
}

clAssign <- function(U){
	apply(U, 1, function(u) order(u, decreasing = TRUE)[1])
}

clusterHeatMap <- function(V){
	k <- ncol(V)
	o <- apply(V, 1, function(v) length(unique(round(v,6))))
	
	vdf <- data.frame(Cluster = rep(1:k, each = nrow(V)), 
			Variable = rep(rownames(V)[order(o)], times = k), 
			Value = c(V[order(o),]))
	vdf$Variable <- factor(vdf$Variable, 
			levels = rownames(V)[order(o)])
	vdf$Cluster <- as.factor(vdf$Cluster)
	p <- ggplot(vdf, aes(Variable, Cluster)) + 
			geom_tile(aes(fill = Value), color = 'white') + 
			scale_fill_gradient(low = 'blue', high = 'yellow')
	p + theme(axis.ticks = element_blank(), 
			axis.text.x = element_text(size = 5, 
					angle = 270, hjust = 0, colour = "grey50"))
	
}

