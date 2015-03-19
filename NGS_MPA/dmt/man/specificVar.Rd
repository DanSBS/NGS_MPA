\name{specificVar}
\alias{specificVar}

\title{Data-specific variation retained in the combined drCCA representation}
\description{
 A function for estimating the amount of data-set specific variation 
 (i.e. variation that is not present in any of the other data sets)
 retained in the combined data set of given dimensionality.
}

\usage{specificVar(datasets,regcca,dim,pca=FALSE)}
\arguments{
  \item{datasets}{A list containing the data matrices to be combined.
        Each matrix needs to have the same number
        of rows (samples), but the number of columns (features) can differ.
        Each row needs to correspond to the same sample in every matrix.}
  \item{regcca}{Output of \code{regCCA} function, containing the
	solution of the generalized CCA.}
  \item{dim}{ The number of dimensions of projected data to be used}
  \item{pca}{A logical variable with default value FALSE. If the value is TRUE,
        the data-specific variation will also be calculated for the PCA
	projected data, where PCA is performed on the columnwise concatenation
	of the given data sets.}
}
\details{
  The function estimates the amount of data-specific information retained in a
  previously calculated drCCA solution. The function uses SVD 
  to estimate the variance of each data set in the drCCA projection of the 
  given dimensions. Data-specific variance is defined as the sum of singular
  values for the covariance matrix of a data set. The value is normalized 
  so that the variation for each of the original data sets is 1. The
  average of the data-specific variances in the projection is also
  calculated. A solution truly focusing on the dependencies usually has 
  a value that grows roughly linearly when the number of dimensions is 
  increased.
  The function can also be used to estimate the same quantity for simple
  PCA projection of the concatenation of the data sets. This can be used as
  a comparison value. For details, please check the reference.
}

\value{
 The function returns a list of following values

   \item{cc}{ Data Specific variation for a drCCA projection of given
              number of dimensions
            }
   \item{pc}{A vector containing the data-specific variations for a PCA
             projection of given dimensions, if pca = TRUE is given
            }
   \item{mcca}{Mean of data-specific variations for a drCCA projection}
   \item{mpca}{Mean of data-specific variation for a PCA projection,
               if pca = TRUE is given}
  

}

\references{ 
Tripathi A., Klami A., Kaski S. (2007), Simple integrative preprocessing preserves what is shared in data sources.
}

\author{Abhishek Tripathi \email{abhishektripathi.at@gmail.com}, Arto Klami}


\seealso{\code{\link{sharedVar}}}

\examples{

#       data(expdata1)
#       data(expdata2)
#       r <- regCCA(list(expdata1,expdata2))

#       specificVar(list(expdata1,expdata2),r,4)


}
\keyword{multivariate}

