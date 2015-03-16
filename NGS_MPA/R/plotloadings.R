# Function to plot CCA projections and vectors
# 
# Author: Daniel V. Samarov
###############################################################################
plotloadings <- function(proj, load, d, ...){
	par(mfrow = c(d, d))
	d2 <- (d^2)
	mat <- matrix(1:d2, d, d, byrow = TRUE)
	mat[lower.tri(mat)] <- NA
	indx <- c(t(mat))
	dindx <- diag(mat)
	
	pltindx <- as.matrix(expand.grid(1:d, 1:d))[,2:1]
	for(i in 1:d2){
		if(is.na(indx[i]))
			plotnull()
		else{
			if(i == dindx)
			ind <- pltindx[i, ]
			
		}
	}
}

plotnull <- function(){
	plot(0, 0, xaxt = 'n', yaxt = 'n',
			xlab = '', ylab = '',
			bty = 'n')
}