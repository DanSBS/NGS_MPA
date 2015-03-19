\name{centerData}
\alias{centerData}
\title{Center data matrix.}
\description{Center data matrix to 0 for each variable by removing the means.}
\usage{
centerData(X, rm.na = TRUE, meanvalue = NULL)
}
\arguments{
  \item{X}{The data set: samples x features. Each feature will be
centered.}
  \item{rm.na}{Ignore NAs.}
  \item{meanvalue}{Can be used to set a desired center value. The
default is 0.}
}
%\details{}
\value{Centered data matrix.}
\references{See citation("dmt").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\note{Note that the model assumes samples x features matrix, and centers each feature.}
%\seealso{}
\examples{
#centerData(X)
}
\keyword{ utilities }
\keyword{ maths}

