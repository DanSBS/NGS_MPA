\name{compare.estimate.and.truth}
\alias{compare.estimate.and.truth}
\title{Measuring model accuracy.}
\description{Compare estimated parameters and the original parameters 
               used to generate the toydata.}
\usage{compare.estimate.and.truth(res, toy)}
\arguments{
  \item{res}{DependencyModel object.}
  \item{toy}{Toydata which was used to learn the model.} 
}
\details{Given original data (including the original parameters) and 
       model learned based on the data, compare the learned parameters to the
       original ones.}
\note{Additional tests added later.}
\value{
Vector with following elements
  \item{wtw.x, wtw.y}{Correlation between the original and estimated values for W*t(W) for X and Y, respectively.}
  \item{phi.x, phi.y}{Correlation between the original and estimated values for phi for X and Y, respectively.}
}
\references{See citation("dmt") for references.}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{\code{\link{generate.toydata}}}
\examples{ 
# toy <- generate.toydata(N = N, 
#           zDim = zdim, xDim = xdim, yDim = ydim,
#           marginal.covariances = marginalCovariances,
#           priors = priors)
# res <- fit.dependency.model(toy$X, toy$Y,
#                      zDimension = zdim,
#              marginalCovariances = marginalCovariances,
#		      priors = priors,
#                     matched = FALSE)
# vec <- compare.estimate.and.truth(res, toy)
#
}
\keyword{ utilities }
