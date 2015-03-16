# Analysis of multi-platform NGS data
# 
# Author: Daniel V. Samarov
###############################################################################
library(PMA)

d <- dir('data')
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
	xtmp <- scale(apply(as.matrix(tmp[,-s]), 2, asinh))
	X[[i]] <- xtmp
}

## Run sparse CCA
cca <- MultiCCA(X, 2, ncomponents = 6, trace = TRUE, standardize = FALSE)
plist <- lapply(1:nf, function(i) X[[i]] %*% cca$ws[[i]])

