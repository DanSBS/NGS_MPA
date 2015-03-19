\name{drCCAcombine}
\alias{drCCAcombine}
\title{A function to combine several data sets}
\description{

Performs drCCA on a collection of data sets with co-occurring samples.
The method utilizes regularized canonical correlation analysis to find
linear projections for each of the data sets, and uses those to
construct a combined representation of lower dimensionality than the
original collection. The method suggests a specific dimensionality for
the combined representation, but it is possible to obtain also
combined data sets of different dimensionality.

}
\usage{drCCAcombine(datasets, reg=0, nfold=3, nrand=50)}

\arguments{
  \item{datasets}{A list containing the data matrices to be combined.
        Each matrix needs to have the same number
	of rows (samples), but the number of columns (features) can differ.
	Each row needs to correspond to the same sample in every matrix.}
  \item{reg}{

Regularization parameter for the whitening step used to remove
data-set specific variation. The value of parameter must be between 0
and 1. The default value is set to 0, which means no regularization
will be used. If a non-zero value is given it means that some of the
dimensions with the lowest variance are ignored when whitening. In
more technical terms, the dimensions whose total contribution to the
sum of eigenvalues of the covariance matrix of each data set is below
reg will not be used for the whitening.

}
  \item{nfold}{The number of cross-validation folds used in the automatic
     dimensionality estimation process. The default value is 3.}
  \item{nrand}{The number of random comparison data-sets created
     for the automatic dimensionality estimation process. The default
     value is 50.}

}
\details{

The function uses \code{\link{regCCA}} to perform the canonical
correlation analysis. The dimensionality of the combined data set is
selected using a statistical test that aims to find which dimensions
capture shared variation significantly more than what would be
found under the assumption that the data sets were independent. For
this purpose rnand collections of random matrices with similar
variance structure but no between-data dependencies are created. The
amount of variation each dimension extracts from leave-out data in the
cross-validation setting with nfold folds is compared to the
distribution obtained from the random matrices, and the dimensions
that differ significantly from the null hypothesis of independence are
kept in the combined representation. For details, please check the 
reference.
 
 }

\value{

  The function returns a list of two values.
 
  \item{proj}{The representation obtained by combining the source data sets.
 This is a matrix that contains a feature representation for each of the
 samples in the analyzed collection. Each row in this result matches the
 corresponding row in the original data sets.}
  \item{n}{The number of dimensions in the combined representation. This is
 equal to ncol(proj).}

}

\references{Tripathi A., Klami A., Kaski S. (2007), Simple integrative
preprocessing preserves what is shared in data sources.}

\author{Abhishek Tripathi \email{abhishektripathi.at@gmail.com}, Arto Klami}

\seealso{\code{\link{regCCA}}}
\examples{

    # data(expdata1)
    # data(expdata2)
    # drCCAcombine(list(expdata1,expdata2),0,2,3)

}

\keyword{multivariate}% at least one, from doc/KEYWORDS
