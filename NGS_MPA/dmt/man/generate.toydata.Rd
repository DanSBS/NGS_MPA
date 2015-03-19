\name{generate.toydata}
\alias{generate.toydata}
\title{Toy data generation}
\description{Generate simulated data which follows the distributional assumptions of the model.}
\usage{generate.toydata(N = 100, zDim = 2, xDim = 3, yDim = 3, 
          marginal.covariances = "full", priors = NULL)}
\arguments{
  \item{N}{ Sample size. }
  \item{zDim}{ Dimensionality of the latent variable. }
  \item{xDim}{ Dimensionality of X data set. }
  \item{yDim}{ Dimensionality of Y data set. }
  \item{marginal.covariances}{"full": full covariance matrices for marginal noise (assumed by pCCA); "diagonal": diagonal covariances for marginal noise (pFA); "isotropic": isotropic covariances (pPCA). }
  \item{priors}{ Set priors for toydata generation. Use as in \code{\link{fit.dependency.model}}.}
}
\details{ Assuming normally distributed latent variables for shared component Z, and data-specific components Zx, Zy. These follow standard multivariate normal distribution N(0, I). The observations X and Y are obtained as X = Wx*Z + Bx*Zx, Y = Wy*Z + By*Zy.}
\value{
List with the following components:
  \item{Z, Zx, Zy }{Shared and data-set specific latent variables.}
  \item{Wx, Wy, Bx, By }{ Transformation matrices. }
  \item{X, Y }{ Data sets. }
}
\references{See citation("dmt") for references.}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\seealso{\code{\link{help}}}
\examples{ toy <- generate.toydata(N = 100, 
                   zDim = 1, xDim = 3, yDim = 3, 
                   marginal.covariances = "full") 
}
\keyword{ utilities }
