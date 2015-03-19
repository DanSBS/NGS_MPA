\name{regCCA}
\alias{regCCA}
\title{Generalized Canonical Correlation Analysis }
\description{
 Solve generalized CCA. Contains a possibility to regularize
 the solution to reduce the effect of noise.
}
\usage{
regCCA(datasets, reg=0)
}
\arguments{
  \item{datasets}{A list containing the data matrices to be analyzed.
        Each matrix needs to have the same number of rows (samples),
        but the number of columns (features) can differ.
	Each row needs to correspond to the same sample in every matrix.}
  \item{reg}{

Regularization parameter for the whitening step used to remove
data-set specific variation. The value of parameter must be between 0
and 1. The default value is set to 0, which means no regularization
will be used. If a non-zero value is given it means that some of the
dimensions with the lowest variance are ignored when whitening. In
more terms, the dimensions whose total contribution to sum of
eigenvalues of the covariance matrix of each data set is below reg
will not be used for the whitening.}  }

\details{
  The function implements generalized CCA by explicitly whitening the
  data sets and then performing a principal component analysis on
  the collection of whitened data sets, instead of directly solving
  the generalized eigenproblem. Singular value decomposition
  is used for both the whitening and the PCA phase, and row-wise mean values of
  each data set are removed before whitening.
}
\value{
  The function returns a list with following 
  components
  
  \item{eigval}{Generalized canonical correlations. In case of two data
	sets (eigval-1) would give the correlations.}
  \item{eigvecs}{List of projection matrices, one for each data set. Each
	projection matrix is a N times m matrix where N is the number
	of samples and m is the total number of dimensions in all of the
	data sets.}
  \item{proj}{Projection of the original data sets by the
	corresponding projection matrices.}

  \item{meanvec}{An array containing columnwise mean vectors for each data matrix}
  \item{white}{An array of whitening matrices for each data set. This might not be of user interest but this value is  used as input in other functions in the package.}

  The function also prints whether regularization was used or not.

}
\references{ 
	Hotelling H. (1936), Relations between two sets of variables, \emph{Biometrika}, \bold{28}, 321-327.

        Kettenring J.R. (1971), Canonical Analysis of several sets of variables, \emph{Biometrika}, \bold{58:3}, 433-451.

       	Tripathi A., Klami A., Kaski S. (2007), Simple integrative
preprocessing preserves what is shared in data sources.
}


\author{Abhishek Tripathi \email{abhishektripathi.at@gmail.com}, Arto Klami}

\seealso{\code{cancor},\code{prcomp},\code{svd} }
\examples{

#     data(expdata1)
#     data(expdata2)

     #performing regCCA
#    test <- regCCA(list(expdata1,expdata2),0) #list of result is stored in test
 
#     test$eigval #generalized canonical correlations
#     test$eigvecs #gCCA components
#     test$proj #projection of data onto gCCA components
#     test$meanvec #array of columnwise mean vectors for each matrix
#     test$white # array of whitening matrix
}

\keyword{multivariate}% at least one, from doc/KEYWORDS

