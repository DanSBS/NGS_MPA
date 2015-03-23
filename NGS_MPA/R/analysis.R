# Analysis of multi-platform NGS data
# 
# Author: Daniel V. Samarov
###############################################################################
library(PMA)

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

## Run sparse CCA
perm <- MultiCCA.permute(X)
clus <- kmeans(r$proj[,1:8], 8)$clus
mm <- model.matrix(~-1+as.factor(group))
X[[5]] <- mm
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

i=1
for(i in 1:nf){
	png(paste('plots/sparse_cca',i,'.png',sep=''), height = 800, width = 800)
	
	plotcca(plist[[i]],cca$ws[[i]],4,
			colnames(X[[i]]),pch=19,cex=0.5, col = group,
			pltname = d[i], highlight = NULL, cex.main = 2)
	
	dev.off()
}
