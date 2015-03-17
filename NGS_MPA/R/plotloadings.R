# Function to plot CCA projections and vectors
# 
# Author: Daniel V. Samarov
###############################################################################
plotcca <- function(proj, load, d, names, 
		pltname = NULL, highlight = NULL, ...){
	par(mfrow = c(d, d))
	par(mar = c(1,2.5,3,1))
	d2 <- (d^2)
	mat <- matrix(1:d2, d, d, byrow = TRUE)
	mat[lower.tri(mat)] <- NA
	indx <- c(t(mat))
	dindx <- diag(mat)
	
	pltindx <- as.matrix(expand.grid(1:d, 1:d))[,2:1]
	
	for(i in 1:d2){
		ind <- pltindx[i, ]
		if(is.na(indx[i]))
			plotnull()
		else{
			if(ind[1] == ind[2]){
				loadi <- load[,ind[1]]
				namesi <- names[loadi != 0]
				loadi <- loadi[loadi != 0]
				aloadi <- abs(loadi)
				cols <- rep('green', length(loadi))
				cols[loadi <= 0] <- 'red'
				o <- order(aloadi, decreasing = TRUE)
				m <- min(5,length(o))
				o <- o[1:m]
				b <- barplot(aloadi[o], 
						las = 2, col = cols[o],
						main = paste('Loading', ind[1]),
						cex.main = 2)
				text(b,max(aloadi)/2,namesi[o],srt=90,cex=1.5)
			}
			else{
				plot(proj[,ind[2]],proj[,ind[1]],
						xlab = '', ylab = '', ...)
				if(!is.null(highlight)){
					cl <- list(...)$col
					pch <- list(...)$pch
					cex <- list(...)$cex
					
					if(!is.null(cl)){
						points(proj[,ind[2]][!highlight],
								proj[,ind[1]][!highlight],
								col = cl[!highlight],
								pch = pch,
								cex = cex)
						points(proj[,ind[2]][highlight],
								proj[,ind[1]][highlight],
								col = cl[highlight],
								pch = pch,
								cex = cex)
					}
				}
			}
		}
		if(!is.null(pltname) & (ind[1] == d) & (ind[2] ==1)){
			text(1,1,label=pltname, cex = 2)
		}
	}
}

plotnull <- function(){
	plot(0:2, 0:2, xaxt = 'n', yaxt = 'n',
			type = 'n',
			xlab = '', ylab = '',
			bty = 'n')
}