\name{sharedVar}
\alias{sharedVar}

\title{Shared variation retained in the combined drCCA representation}
\description{
 A function for estimating the amount of shared variation (i.e. variation
 that is common to more than one data set) retained in the combined data 
 set of given dimensionality.
}
\usage{sharedVar(datasets,regcca,dimension,pca=FALSE)}
\arguments{
  \item{datasets}{A list containing the data matrices to be combined.
        Each matrix needs to have the same number
        of rows (samples), but the number of columns (features) can differ.
        Each row needs to correspond to the same sample in every matrix.}
  \item{regcca}{Output of \code{regCCA} function, containing the
	solution of the generalized CCA.}
  \item{dimension}{The number of dimensions of projected data to be used}
  \item{pca}{A logical variable with default value FALSE. If the value is TRUE,
        the pairwise variation will also be calculated for the PCA projected 
        data, where PCA is performed on the columnwise concatenation of the
	given data sets.}
}

\details{
  The function estimates the amount of shared information retained in a
  previously calculated drCCA solution. It calculates the shared
  variation between all pairs of the data sets returned from drCCA
  combined data for a particular dimensionality. The function also
  calculates the same quantities for the original data and for the simple
  PCA projection of the concatenation of data sets. This can be used as a 
  comparison value. If the full dimensionality of drCCA projection or
  the PCA projection is used, the sum of all pairs of shared variations will
  be the same. The mean of shared variations for drCCA and PCA 
  is estimated, normalized in a way that the value for original data
  sets will be 1. A good result will have value greater than 1. For
  details please refer to the reference below.
 

}

\value{
  A list of following elements is returned

    \item{oo}{A matrix containing the pairwise shared variations for original 
             data sets}
   \item{cc}{A matrix containing the pairwise shared variations for a drCCA
             projection of given dimensions
            }
   \item{pc}{A matrix containing the pairwise shared variations for a PCA
             projection of given dimensions, if pca = TRUE is given
            }
   \item{mcca}{Mean of shared variation between all pairs for drCCA}
   \item{mpca}{Mean of shared variation between all pairs for PCA,
               if pca = TRUE is given}
}

\references{ 
	Tripathi A., Klami A., Kaski S. (2007), Simple integrative
preprocessing preserves what is shared in data sources.

}
\author{Abhishek Tripathi \email{abhishektripathi.at@gmail.com}, Arto Klami}

\seealso{\code{\link{specificVar}}}
\examples{

  #     data(expdata1)
  #     data(expdata2)
  #     r <- regCCA(list(expdata1,expdata2))

  #     sharedVar(list(expdata1,expdata2),r,4)


}
\keyword{multivariate}
