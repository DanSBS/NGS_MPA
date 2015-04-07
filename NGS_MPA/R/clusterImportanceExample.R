# Testing out some concepts
# 
# Daniel V. Samarov
###############################################################################
require(genlasso)
require(MASS)
source("C:/Users/samarov/git/gsmethods/gsmethods/R/smtlasso_v2.R")

c1 <- mvrnorm(200, c(-1,0), matrix(c(.1,.25,.25,.9),2,2))
c2 <- mvrnorm(200, c(1, 0), matrix(c(.1,.25,.25,.9),2,2))
c3 <- mvrnorm(200, c(-1,-3), matrix(c(.1,.25,.25,.9),2,2))
x <- rbind(c1,c2,c3)


k <- 20
ngam <- 25
max.iter <- 100
eps <- 1e-3
lambda <- .1

mfc <- mfClus(x, k, lambda, reduce = TRUE)

mfClus <- function(x, k, 
		lambda, ngam = 100, max.iter = 1000,
		eps = 1e-3, iter = 3,
		reduce = FALSE, rsize = 2000){
	
	## Get initial starting values as centers 
	## from kmeans
	V <- t(kmeans(x,k)$cent)
	
	if(reduce)
		Xred <- kmeans(x, rsize)$cen
	else
		Xred <- x
	
	for(i in 1:iter){
		print(paste('Iteration:', k))
		U <- t(smtlasso(V, t(x), NULL, lambda = lambda, 
						rho = 0, ngam = ngam,
						max.iter = max.iter,
						parallel = FALSE,
						eps = eps, pos = 1, 
						method='sgl', 
						intercept=FALSE,
						singled=1))
		
		
		## Remove 0 columns
		chk0 <- colSums(U) != 0
		U <- U[,chk0]
		
		k <- ncol(U)
		print(paste('New # clusters:', k))
		D <- diag(k) %x% rep(1, (k - 1))
		oD <- c(rep(1, k) %x% (diag((k - 1))*-1))
		D[D != 1] <- oD
		
		if(reduce)
			Ured <- kmeans(U, rsize)$cent
		else
			Ured <- U
		
		Vlist <- lapply(1:ncol(x), function(i){
					gl <- genlasso(Xred[,i], Ured, D,
							btol = 1e-3,
							maxsteps = 500,
							approx = FALSE,
							verbose = TRUE)
					gl
				})
		
		V <- do.call('rbind', lapply(1:length(Vlist), function(i){
							gl <- Vlist[[i]]
							V <- gl$beta
							df <- gl$df
							xi <- x[,i]
							score <- colSums((U %*% V - xi)^2) / df +
									2 * asinh(nrow(x) - df)
							o <- order(score)[1]
							V[,o]
						}))
		
		
	}
	
	rownames(V) <- colnames(x)
	
	return(list(U = U, V = V, k = k))
	
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

pdf('plots/clusterImportance_ILL250_ROC.pdf',
		height = 8, width = 11)
clusterHeatMap(mfc$V)
dev.off()
