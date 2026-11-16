\name{PCAaugHD}
\alias{PCAaugHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-Dimensional Augmentation Estimate for PCA
}
\description{
Estimates the number of spikes in a spiked covariance model using the high-dimensional predictor augmentation (HDPA) estimator of Radojicic and Virta (2025).
In contrast to \code{\link{PCAaug}} the estimate stays consistent when the dimension and the number of augmented components grow proportionally with the sample size.
}
\usage{
PCAaugHD(X, naug = NULL, nrep = 1, sigma2 = NULL, K = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix with \code{n} rows and \code{p} columns. \code{p} can be larger than \code{n}.}
  \item{naug}{number of components to be augmented. If \code{NULL}, \code{max(p, min(ceiling(n/2), 1000))} is used. See details.}
  \item{nrep}{number of repetitions for the augmentation procedure.}
  \item{sigma2}{the noise variance. If \code{NULL} it is estimated from the data, see details.}
  \item{K}{the largest value of k considered. If \code{NULL}, \code{min(20, min(p, n-1) - 1)} is used.}
}
\details{
The model assumes that the eigenvalues of the covariance matrix are of the form \eqn{\lambda_1 \geq ... \geq \lambda_{k} > \lambda_{k+1} =  ... = \lambda_p = \sigma^2}{lambda_1 >= ... >= lambda_k > lambda_k+1 =  ... = \lambda_p = sigma^2}
and the goal is to estimate the value of k, where \eqn{1 \leq k \leq K}{1 <= k <= K}.

Let \eqn{\gamma = p/n}{gamma = p/n} and \eqn{\rho = r/n}{rho = r/n}, where \eqn{r}{r} is the number of augmented components. The augmented estimator of Luo and Li (2021) (see \code{\link{PCAaug}}) is not consistent when
\eqn{\gamma}{gamma} or \eqn{\rho}{rho} is not close to zero. Radojicic and Virta (2025) modify it as follows:
\itemize{
\item The eigenvalues \eqn{l_j}{l_j} of the covariance matrix of the original data are debiased using
\eqn{f(l) = \{l - \sigma^2(1+\gamma)\}/2 + [\{l - \sigma^2(1+\gamma)\}^2 - 4 \sigma^4 \gamma]_+^{1/2}/2}{f(l) = {l - sigma^2(1+gamma)}/2 + [{l - sigma^2(1+gamma)}^2 - 4 sigma^4 gamma]_+^(1/2)/2}.
\item \eqn{r}{r} Gaussian components with variance \eqn{\sigma^2}{sigma^2} are added to the data and \eqn{\|u_j\|^2}{||u_j||^2} denotes the squared norm of the augmented part of the jth eigenvector of the covariance matrix of the augmented data.
\item The rescaled norms \eqn{h_j = f(l_j)\{f(l_j) + (\gamma + \rho)\sigma^2\}/\{f(l_j) + \sigma^2\} \|u_j\|^2}{h_j = f(l_j){f(l_j) + (gamma + rho) sigma^2}/{f(l_j) + sigma^2} ||u_j||^2} are asymptotically constant both for \eqn{j \leq k}{j <= k} and for \eqn{j > k}{j > k}, and the jump between the two levels is negative.
\item The estimate of k is the value \eqn{j = 1, ..., K}{j = 1, ..., K} which minimizes \eqn{h_{j+1} - h_j}{h_{j+1} - h_j}.
}

If \code{sigma2 = NULL}, the noise variance is estimated by the median of the non-zero eigenvalues of the covariance matrix, divided by the median of the Marchenko-Pastur distribution with ratio \eqn{\gamma}{gamma}
if \eqn{\gamma \leq 1}{gamma <= 1}, and by \eqn{\gamma}{gamma} times the median of the Marchenko-Pastur distribution with ratio \eqn{1/\gamma}{1/gamma} otherwise. This estimate is biased upwards when k is not small compared to \code{min(p, n)/2}.
In that case it is better to provide \code{sigma2}.

The jump at k grows with \eqn{\rho}{rho}, so a large number of augmented components is preferable. The computing time grows however with \code{min(n, p + naug)^3}, which is the reason for the cap in the default value of \code{naug}.
If \code{nrep > 1}, the norms \eqn{\|u_j\|^2}{||u_j||^2} are averaged over the repetitions.

The theory assumes Gaussian data and that the spikes satisfy \eqn{\lambda_k - \sigma^2 > \sigma^2(\gamma + \rho)^{1/2}}{lambda_k - sigma^2 > sigma^2 (gamma + rho)^(1/2)}. Radojicic and Virta (2025) note that a sample size of at least 500 is needed for reliable estimates.

The criterion can be visualized using \code{\link{ladleplot}}. Here the index starts at 1 and the element \code{fn} contains \eqn{h_j}{h_j} and \code{phin} the debiased eigenvalues.
}

\value{
A list of class ladle containing:
  \item{method}{the string HDPA.}
  \item{k}{the estimated value of k.}
  \item{fn}{the rescaled norms \eqn{h_j}{h_j}, \eqn{j = 1, ..., K+1}{j = 1, ..., K+1}.}
  \item{phin}{the debiased eigenvalues \eqn{f(l_j)}{f(l_j)}, \eqn{j = 1, ..., K+1}{j = 1, ..., K+1}.}
  \item{gn}{the main criterion, the differences \eqn{h_{j+1} - h_j}{h_{j+1} - h_j}, \eqn{j = 1, ..., K}{j = 1, ..., K}. k is the value where gn takes its minimum.}
  \item{lambda}{the eigenvalues of the covariance matrix.}
  \item{W}{the transformation matrix to the principal components.}
  \item{S}{data matrix with the centered principal components.}
  \item{MU}{the location of the data which was substracted before calculating the principal components.}
  \item{data.name}{the name of the data for which the augmented estimate was computed.}
  \item{sigma2}{the value used as noise variance.}
  \item{index.start}{the index of the first element of the criteria, used by \code{\link{ladleplot}}.}
}
\references{
\cite{Luo, W. and Li, B. (2021), On Order Determination by Predictor Augmentation, Biometrika, 108, 557--574. <doi:10.1093/biomet/asaa077>}

\cite{Radojicic, U. and Virta, J. (2025), Dimension Estimation in a Spiked Covariance Model Using High-Dimensional Data Augmentation, Biometrika, 112, asaf052. <doi:10.1093/biomet/asaf052>}
}
\author{
Klaus Nordhausen
}
\seealso{
\code{\link{PCAaug}}, \code{\link{ladleplot}}, \code{\link{PCAladle}}
}
\examples{
n <- 500
p <- 100
ev <- c(4, 3, 2, rep(0, p - 3)) + 1
X <- matrix(rnorm(n * p), ncol = p) \%*\% diag(sqrt(ev))

test <- PCAaugHD(X)
test
summary(test)
ladleplot(test)
ladleplot(test, crit = "fn")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
