\name{plotVar}
\alias{plotVar}

\title{Data-specific and shared variance for several dimensionalities}
\description{
  A function for calculating the captured variations for several
  different number of retained dimensions. This is a wrapper
  over \code{specificVar} and \code{sharedVar} to help
  computing both for a range of dimensionalities.
}
\usage{plotVar(datasets,regcca,dimVector,plot=FALSE)}
\arguments{
  \item{datasets}{A list containing the data matrices to be combined.
        Each matrix needs to have the same number
        of rows (samples), but the number of columns (features) can differ.
        Each row needs to correspond to the same sample in every matrix.}
  \item{regcca}{Output of \code{regCCA} function, containing the
	solution of the generalized CCA.}
  \item{dimVector}{A list of dimensions for which the retained variations
	are to be computed.}
  \item{plot}{A logical variable with default value FALSE. If the value is
	 TRUE, the functions creates a plot of the output.}
}
\details{
The function uses \code{specificVar} and \code{sharedVar}
to do all the computation. The purpose of this function is to provide an
easy way to visualize the properties of the reduced-dimensional
representation created by drCCA. The function also estimates the
same quantities for PCA of concatenated feature vectors to illustrate
the difference to optimal linear model based on preserving the
total variation in the whole collection of data sets.
}

\value{ 
        The function returns the data-specific and shared variance
  for the given values of dimensions in a list. The list has four 
  components.
  
  \item{pw_cca}{A vector with values as shared variances captured by drCCA
                for the given dimensions}

  \item{pw_pca}{A vector with values as shared variances captured by PCA
                for the given dimensions}

  \item{within_cca}{A vector with values as data-specific variances captured 
                    by drCCA for the given dimensions}

  \item{within_pca}{A vector with values as data-specific variances captured 
                    by PCA for the given dimensions}
}

\references{ 
	Tripathi A., Klami A., Kaski S. (2007), Simple integrative preprocessing preserves what is shared in data sources, \emph{submitted for publication}.

}

\author{Abhishek Tripathi, Arto Klami}


\seealso{\code{\link{sharedVar}},\code{\link{specificVar}}}

\examples{

#       data(expdata1)
#       data(expdata2)
#       r <- regCCA(list(expdata1,expdata2))#

#       plotVar(list(expdata1,expdata2),r,c(1:2),4)


}
\keyword{multivariate}

