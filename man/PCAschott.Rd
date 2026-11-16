\name{PCAschott}
\alias{PCAschott}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing for Subsphericity using Schott's High-Dimensional Test
}
\description{
Tests the equality of the last p-k eigenvalues of the covariance matrix using the test of Schott (2006), which is valid also when the dimension is large compared to the sample size.
}
\usage{
PCAschott(X, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix with p>1 columns. The number of columns can be larger than the number of rows.}
  \item{k}{the number of eigenvalues larger than the equal ones. Can be between 0 and p-2.}
}
\details{
Let \eqn{l_1 \geq \ldots \geq l_p}{l_1 >= ... >= l_p} be the eigenvalues of the sample covariance matrix of the n observations and \eqn{q = p-k}{q = p-k}.
The test statistic is based on the generalization of the sphericity statistic of Ledoit and Wolf (2002)
\deqn{U = \frac{q^{-1}\sum_{i=k+1}^p l_i^2}{(q^{-1}\sum_{i=k+1}^p l_i)^2} - 1}{U = (sum_{i=k+1}^p l_i^2 / q) / (sum_{i=k+1}^p l_i / q)^2 - 1}
and is \eqn{T = (n-1-k)qU/2}{T = (n-1-k) q U / 2}. Under the null hypothesis and for normal data, \eqn{T} has asymptotically a chisquare distribution with \eqn{q(q+1)/2-1}{q(q+1)/2-1} degrees of freedom.
Schott (2006) shows that this approximation holds both when n goes to infinity for a fixed dimension and when n and p go to infinity together, also for p larger than n.

For data where p is small compared to n, also the likelihood ratio test \code{\link{PCAlrt}} can be used. Its chisquare approximation becomes however poor when p is not small compared to n and it requires n > p.
}
\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the degrees of freedom of the test.}
  \item{method}{character string which test was performed.}
  \item{data.name}{character string giving the name of the data.}
  \item{alternative}{character string specifying the alternative hypothesis.}
  \item{k}{the number or larger eigenvalues used in the testing problem.}
  \item{W}{the transformation matrix to the principal components.}
  \item{S}{data matrix with the centered principal components.}
  \item{D}{the underlying eigenvalues.}
  \item{MU}{the mean vector of the data which was substracted before calculating the principal components.}
  \item{SCATTER}{the computed covariance matrix matrix.}
}
\references{
\cite{Schott, J.R. (2006), A High-Dimensional Test for the Equality of the Smallest Eigenvalues of a Covariance Matrix, Journal of Multivariate Analysis, 97, 827--843. <doi:10.1016/j.jmva.2005.05.003>}

\cite{Ledoit, O. and Wolf, M. (2002), Some Hypothesis Tests for the Covariance Matrix when the Dimension is Large Compared to the Sample Size, The Annals of Statistics, 30, 1081--1102. <doi:10.1214/aos/1031689018>}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{PCAlrt}}, \code{\link{PCAasymp}}, \code{\link{PCAboot}}, \code{\link{PCAaugHD}}
}
\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))
PCAschott(X, 2)

# more variables than observations
n <- 50
p <- 100
X2 <- matrix(rnorm(n * p), ncol = p) \%*\% diag(c(5, 4, rep(1, p - 2)))
PCAschott(X2, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }
