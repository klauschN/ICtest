\name{PCAlrt}
\alias{PCAlrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing for Subsphericity using the Likelihood Ratio Test
}
\description{
Tests the equality of the last p-k eigenvalues of the covariance matrix using the likelihood ratio test for normal data with the Bartlett-Lawley correction.
}
\usage{
PCAlrt(X, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix with p>1 columns and more rows than columns.}
  \item{k}{the number of eigenvalues larger than the equal ones. Can be between 0 and p-2.}
}
\details{
Let \eqn{l_1 \geq \ldots \geq l_p}{l_1 >= ... >= l_p} be the eigenvalues of the sample covariance matrix of the n observations, \eqn{q = p-k}{q = p-k} and \eqn{\bar l}{l.bar} the mean of the last q eigenvalues.
The test statistic is
\deqn{T = c_n \left(q \log \bar l - \sum_{i=k+1}^p \log l_i\right)}{T = c_n (q log(l.bar) - sum_{i=k+1}^p log(l_i))}
with the Bartlett-Lawley correction factor
\deqn{c_n = n - 1 - k - \frac{1}{6}\left(2q + 1 + \frac{2}{q}\right) + \bar l^2 \sum_{i=1}^k \frac{1}{(l_i - \bar l)^2}.}{c_n = n - 1 - k - (2q + 1 + 2/q)/6 + l.bar^2 sum_{i=1}^k 1/(l_i - l.bar)^2.}
Under the null hypothesis and for normal data, \eqn{T} has asymptotically a chisquare distribution with \eqn{q(q+1)/2-1}{q(q+1)/2-1} degrees of freedom when n goes to infinity.

The test requires n > p. The chisquare approximation is good when p is small compared to n but the test becomes increasingly liberal when p grows with n, see Schott (2006).
In that case \code{\link{PCAschott}} should be used.
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
\cite{Bartlett, M.S. (1950), Tests of Significance in Factor Analysis, British Journal of Statistical Psychology, 3, 77--85. <doi:10.1111/j.2044-8317.1950.tb00285.x>}

\cite{Lawley, D.N. (1956), Tests of Significance for the Latent Roots of Covariance and Correlation Matrices, Biometrika, 43, 128--136. <doi:10.1093/biomet/43.1-2.128>}

\cite{Schott, J.R. (2006), A High-Dimensional Test for the Equality of the Smallest Eigenvalues of a Covariance Matrix, Journal of Multivariate Analysis, 97, 827--843. <doi:10.1016/j.jmva.2005.05.003>}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{PCAschott}}, \code{\link{PCAasymp}}, \code{\link{PCAboot}}
}
\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))
PCAlrt(X, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }
