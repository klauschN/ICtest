\name{SIRaug}
\alias{SIRaug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Augmentation Estimate for SIR
}
\description{
In the supervised dimension reduction context with response y and explaining variables x, this function provides the augmentation estimate
for the dimension of the central subspace for SIR.
}
\usage{
SIRaug(X, y, h = 10, naug = 1, nrep = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix.}
  \item{y}{numeric response vector.}
  \item{h}{number of slices in SIR.}
  \item{naug}{number of components to be augmented.}
  \item{nrep}{number of repetitions for the augmentation procedure.}
  \item{\dots}{arguments passed on to \code{\link{quantile}}.}
}
\details{
The idea here is that the eigenvalues of the SIR-M matrix are of the form  \eqn{\lambda_1 \geq ... \geq \lambda_k > 0 =  ... = 0}{lambda_1 >= ... >= lambda_k > 0 =  ... = 0} and the eigenvectors
of the non-zero eigenvalue span the central subspace.

The augmented estimator adds for that purpose \code{naug} standard Gaussian components which are independent of the response. As SIR is affine equivariant the variance of the
augmented components does not matter and no estimation of a noise variance is needed. The augmented components belong then to the part of the predictors which is
not related to the response and the augmentation estimator uses them to measure the variation of the eigenvectors. For a more stable result it is recommended to repeat
the augmentation process several times and Luo and Li (2021) recommend to use for \code{naug} approximately \code{p/5} or \code{p/10} where \code{p} is the number of columns of \code{X}.

The augmented estimator combines then the values of the scaled eigenvalues and the variation measured via augmentation. The main idea there is that for
distinct eigenvalues the variation of the eigenvectors is small and for equal eigenvalues the corresponding eigenvectors have large variation.

The augmented estimate for k is the value where the measure takes its minimum and can be also visualized as a ladle.

For further details see Luo and Li (2021).
}

\value{
A list of class ladle containing:
  \item{method}{the string SIR.}
  \item{k}{the estimated value of k.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the augmented components for the different number of components.}
  \item{phin}{normalized eigenvalues of the M matrix.}
  \item{gn}{the main criterion for the augmented estimate - the sum of fn and phin. k is the value where gn takes its minimum}
  \item{lambda}{the eigenvalues of the M matrix.}
  \item{W}{the transformation matrix to the supervised components.}
  \item{S}{data matrix with the centered supervised components.}
  \item{MU}{the location of the data which was substracted before calculating the supervised components.}
  \item{data.name}{the name of the data for which the augmented estimate was computed.}
}
\references{
\cite{Luo, W. and Li, B. (2021), On Order Determination by Predictor Augmentation, Biometrika, 108, 557--574. <doi:10.1093/biomet/asaa077>}
}
\author{
Klaus Nordhausen
}
\seealso{
\code{\link{ladleplot}}, \code{\link{SIRladle}}, \code{\link{PCAaug}}, \code{\link{FOBIaug}}
}
\examples{
n <- 1000
p <- 10
A <- matrix(rnorm(p*p), ncol=p)
x <- matrix(rnorm(n*p), ncol=p)
eps <- rnorm(n, sd=0.1)
y <- x[,1] / (0.5 + (x[,2] + 1.5)^2) + eps
X <- x \%*\% t(A)

test <- SIRaug(X, y, nrep = 10)
test
summary(test)
plot(test)
ladleplot(test)
ladleplot(test, crit = "fn")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
