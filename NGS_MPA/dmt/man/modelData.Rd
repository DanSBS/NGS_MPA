\name{modelData}
\alias{X}
\alias{Y}
\docType{data}
\title{ An example data set }
\description{
            Preprocessed gene expression and gene copy number levels of 51 patients in 
            chromosome 17 for 10 genes.
}
\usage{data(modelData)}
\format{
  \describe{
        \item{X}{Gene expression levels in matrix form. Genes are in
      columns and samples in rows.}
      
      \item{Y}{Gene copy number levels in matrix form. Genes are in
      columns and samples in rows.}
  }
}
\source{ Integrated gene copy number and expression microarray analysis
  of gastric cancer highlights potential target genes.  Myllykangas et
  al., \emph{International Journal of Cancer}, vol. \bold{123},
  \bold{no. 4}, pp. 817--25, 2008.}
\examples{
data(modelData)

model <- fit.dependency.model(X, Y)
model

}
\keyword{datasets}
