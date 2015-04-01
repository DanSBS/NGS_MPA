# Analysis of multi-platform NGS data
# 
# Author: Daniel V. Samarov
###############################################################################
library(PMA)
source('R/plotcca.R')
source('R/mvcca.R')

d <- dir('data')
d <- d[d!='PlatGen_ROC_wclus.csv']
nf <- length(d)

X <- vector('list',nf)

for(i in 1:nf){
	f <- paste('data/',d[i],sep='')
	print(f)
	tmp <- read.csv(f,as.is=TRUE,
			stringsAsFactors=FALSE)
	
	if(d[i] == "PlatGen_ROC.csv")
		tmp <- tmp[,names(tmp) != 'R_Ins_10_percentile']
	s <- 1:4
	nc <- ncol(tmp)
	s <- c(s, (nc-10):nc)
	xtmp <- scale(apply(as.matrix(tmp[tmp$GROUP %in% c(1:7),-s]), 2, asinh))
	xtmp <- xtmp[,!apply(xtmp,2,function(u) all(is.nan(u)))]
	X[[i]] <- xtmp
}

group <- read.csv('data/PlatGen_ROC_wclus.csv',as.is=TRUE,
		stringsAsFactors=FALSE)$Cluster

cCCA <- clusCCA(X, group, 0.1, frac = 0.1)
cCCA.pen <- clusCCA.pen(X, group, 1)
ncols <- unlist(lapply(X, ncol))


## Run sparse CCA

cca <- MultiCCA(X, 1.5, ncomponents = 4, 
		trace = TRUE, standardize = FALSE)
plist <- lapply(1:nf, function(i) X[[i]] %*% cca$ws[[i]])

## Generating plots highlighting groups 1, 6, and 13 as well
## as variables selected in first component
cols <- rep(rgb(190,190,190,50,maxColorValue=255), nrow(tmp))
cols[tmp$GROUP == 1] <- 'blue'
cols[tmp$GROUP == 6] <- 'green'
cols[tmp$GROUP == 13] <- 'red'

highlight <- rep(FALSE, nrow(tmp))
highlight[tmp$GROUP %in% c(1,6,13)] <- TRUE

i=2
for(i in 1:nf){
	png(paste('plots/sparse_cca',i,'.png',sep=''), height = 800, width = 800)
	
	plotcca(p[,1:ncol(X[[i]])],cCCA[1:ncol(X[[i]]),1:ncol(X[[i]])],4,
			colnames(X[[i]]),pch=group,cex=0.5, col = group,
			pltname = d[i], highlight = NULL, cex.main = 2)
	
	dev.off()
}

k <- 1
par(mfrow=c(2,2))
for(i in 1:4)
{
	if(i != 1)
		k <- k + ncols[i-1]
	pPen <- X[[i]] %*% cCCA.pen[k:sum(ncols[1:i]), 1:4]
	
	png(paste('plots/cluster_cca1_v2_', i, '.png', sep =''), height = 800, width = 800)
	plotcca(pPen[,1:4],cCCA.pen[k:sum(ncols[1:i]), 1:4],4,
			colnames(X[[i]]),pch=group,cex=0.5, col = group,
			pltname = d[i], highlight = NULL, cex.main = 2)
	dev.off()
}
