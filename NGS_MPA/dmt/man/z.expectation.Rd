\name{z.expectation}
\alias{z.expectation}
\title{ Expectation of the latent variable. }
\description{ Calculates expectation of the latent variable, given
data and model parameters.}
\usage{
z.expectation(model, X, Y = NULL)
}
\arguments{
  \item{model}{ \linkS4class{DependencyModel} }
  \item{X}{ Data set X. }
  \item{Y}{ Optional, second data set. Used in dependency models where
two co-occurring observations are assumed to stem from a shared latent
variable. }
}
%\details{}
\value{ Returns a matrix (latent features x samples). }
\references{ See citation("dmt") for references. }
\author{ Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%\section{Warning }{....}
%\seealso{\code{\link{help}}}
\examples{
library(dmt)
data(modelData)
res <- pfa(X, Y, zDimension = 2)
}
\keyword{ utilities }

