\name{fit.dependency.model}
%\Rdversion{1.1}
\alias{fit.dependency.model}
\alias{ppca}
\alias{pcca}
\alias{pfa}
%\alias{pcca.isotropic}
\title{Fit dependency model between two data sets.}

\description{Fit generative latent variable model (see vignette for
model specification) on two data sets. Regularize the solutions with
priors, including constraints on marginal covariance structures, the
structure of W, latent dimensionality etc. Probabilistic versions of
PCA, factor analysis and CCA are available as special cases.}

\usage{
fit.dependency.model(X, Y, zDimension = 1, marginalCovariances = "full",
                     epsilon = 1e-3,
                     priors = list(), matched = TRUE,
                     includeData = TRUE, calculateZ = TRUE, verbose = FALSE)
ppca(X, Y = NULL, zDimension = NULL, includeData = TRUE, calculateZ = TRUE)
pfa(X, Y = NULL, zDimension = NULL, includeData = TRUE, calculateZ = TRUE, priors = NULL)
pcca(X, Y, zDimension = NULL, includeData = TRUE, calculateZ = TRUE)
%pcca.isotropic(X, Y, zDimension = NULL, matched = FALSE, epsilon = 1e-6, includeData = TRUE, calculateZ = TRUE)
}
\arguments{
  \item{X, Y}{Data set/s X and Y. 'Variables x samples'. The second data set (\code{Y}) is optional.}
  \item{zDimension}{Dimensionality of the shared latent variable.}
  \item{marginalCovariances}{ Structure of marginal covariances,
assuming multivariate Gaussian distributions for the dataset-specific
effects. Options: \code{"identical isotropic"}, \code{"isotropic"},
\code{"diagonal"} and \code{"full"}. The difference between isotropic
and identical isotropic options is that in isotropic model, phi$X !=
phi$Y in general, whereas with isotropic model phi$X = phi$Y.}
  \item{epsilon}{ Convergence limit. }
  \item{priors}{Prior parameters for the model. A list, which can contain some of the following elements:
    \describe{
  \item{W}{Rate parameter for exponential distribution (should be
positive). Used to specify the prior for Wx and Wy in the dependency
model. The exponential prior is used to produce non-negative solutions
for W; small values of the rate parameter correspond to an
uninformative prior distribution.}
  \item{Nm.wxwy.mean}{ Mean of the matrix normal prior distribution for the
    transformation matrix T. Must be a matrix of size (variables in
    first data set) x (variables in second data set). If value is
    \code{1}, \code{Nm.wxwy.mean} will be made identity matrix of appropriate size.}
  \item{Nm.wxwy.sigma}{ Variance parameter for the matrix normal prior
    distribution of the transformation matrix \code{T}. Described the allowed
    deviation scale of the transformation matrix \code{T} from the mean matrix
    \code{Nm.wxwy.mean}.}
    }}

  \item{matched}{Logical indicating if the variables (dimensions) are
matched between X and Y. Applicable only when dimX = dimY. Affects the
results only when prior on the relationship Wx ~ Wy is set, i.e. when
priors$Nm.wx.wy.sigma < Inf.}
  
  \item{includeData}{Logical indicating whether the original data is
included to the model output. Using \code{FALSE} can be used to
    save memory.}
  
  \item{calculateZ}{Logical indicating whether an expectation of the
latent variable Z is included in the model output. Otherwise the
expectation can be calculated with \code{getZ} or
\code{z.expectation}. Using \code{FALSE} speeds up the calculation of
the dependency model.}
  
  \item{verbose}{ Follow procedure by intermediate messages. }

}

\details{The \code{fit.dependency.model} function fits the dependency
model X = N(W$X * Z, phi$X); Y = N(W$Y * Z, phi$Y) with the
possibility to tune the model structure and parameter priors.

In particular, the dataset-specific covariance structure phi can be
defined; non-negative priors for W are possible; the relation between
W$X and W$Y can be tuned. For a comprehensive set of examples, see the
example scripts in the tests/ directory of this package.

Special cases of the model, obtained with particular prior
assumptions, include probabilistic canonical correlation analysis
(\code{pcca}; \cite{Bach & Jordan 2005}), probabilistic principal
component analysis (\code{ppca}; \cite{Tipping & Bishop 1999}),
probabilistic factor analysis (\code{pfa}; \cite{Rubin & Thayer
1982}), and a regularized version of canonical correlation analysis
(pSimCCA; \cite{Lahti et al. 2009}).
 
The standard probabilistic PCA and factor analysis are methods for a
single data set (X ~ N(WZ, phi)), with isotropic and diagonal
covariance (phi) for pPCA and pFA, respectively. Analogous models for
two data sets are obtained by concatenating the two data sets, and
performing pPCA or pFA.

Such special cases are obtained with the following choices in the
\code{fit.dependency.model} function:

\describe{

\item{pPCA}{ \code{marginalCovariances = "identical isotropic"} (\cite{Tipping & Bishop 1999}) }

\item{pFA}{ \code{marginalCovariances = "diagonal"} (\cite{Rubin & Thayer 1982}) }

\item{pCCA}{ \code{marginalCovariances = "full"} (\cite{Bach & Jordan 2005}) }

\item{pSimCCA}{ \code{marginaCovariances = "full", priors =
list(Nm.wxwy.mean = I, Nm.wxwy.sigma = 0)}. This is the default
method, corresponds to the case with W$X = W$Y.  (\cite{Lahti et al.
2009}) }

\item{pSimCCA with T prior}{ \code{marginalCovariances = "isotropic",
priors = list(Nm.wxwy.mean = 1, Nm.wx.wy.sigma = 1} (\cite{Lahti et al. 2009})
}}

To avoid computational singularities, the covariance matrix phi is
regularised by adding a small constant to the diagonal.

}

\value{ \linkS4class{DependencyModel} }


\references{
Dependency Detection with Similarity Constraints, Lahti et al., 2009
Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal
Processing, \url{http://arxiv.org/abs/1101.5919}

A Probabilistic Interpretation of Canonical Correlation Analysis, Bach
Francis R. and Jordan Michael I. 2005 Technical Report 688. Department
of Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis, Tipping Michael E. and
Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical
Society}, Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}

EM Algorithms for ML Factorial Analysis, Rubin D. and Thayer
D. 1982. \emph{Psychometrika}, \bold{vol. 47}, no. 1.

}
\author{ Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti
\email{leo.lahti@iki.fi} }

\seealso{ Output class for this function:
\linkS4class{DependencyModel}. Special cases: \code{ppca}, \code{pfa}, \code{pcca}}


\examples{
data(modelData) # Load example data X, Y

# probabilistic CCA
model <- pcca(X, Y)

# dependency model with priors (W>=0; Wx = Wy; full marginal covariances)
model <- fit.dependency.model(X, Y, zDimension = 1, 
      	 		      priors = list(W = 1e-3, Nm.wx.wy.sigma = 0), 
			      marginalCovariances = "full")

# Getting the latent variable Z when it has been calculated with the model
#getZ(model)

}
\keyword{math}

