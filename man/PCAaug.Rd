\name{PCAaug}
\alias{PCAaug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Augmentation Estimate for PCA
}
\description{
For p-variate data, the augmentation estimate for PCA assumes that the last p-k eigenvalues are equal. Combining information from the eigenvalues and eigenvectors
of the covariance matrix the augmentation estimator yields an estimate for k. 
}
\usage{
PCAaug(X, noise = "median", naug = 1, nrep = 1, sigma2 = NULL, alpha = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix.}
  \item{noise}{name of the method to be used to estimate the noise variance. Options are \code{"median"}, \code{"last"}, \code{"quantile"} or \code{"known"}. See details.}
  \item{naug}{number of components to be augmented.}
  \item{nrep}{number of repetitions for the augmentation procedure.}
  \item{sigma2}{value of the noise variance when \code{noise = "known"}.}
  \item{alpha}{the quantile to be used when \code{noise = "quantile"}.}
}
\details{
The model here assumes that the eigenvalues of the covariance matrix are of the form \eqn{\lambda_1 \geq ... \geq \lambda_{k} > \lambda_{k+1} =  ... = \lambda_p}{lambda_1 >= ... >= lambda_k > lambda_k+1 =  ... = \lambda_p}
and the goal is to estimate the value of k. The value \eqn{\lambda_{k+1}}{lambda_k+1} corresponds then to the noise variance. 

The augmented estimator adds for that purpose \code{naug} Gaussian components with the provided noise variance which needs to be provided
(\code{noise = "known"}) or estimated from the data. Three estimation methods are available. In the case of \code{noise = "median"} the estimate is the median of the eigenvalues of the covariance matrix, in the case of \code{noise = "last"} it corresponds to the last eigenvalue of the covariance matrix and in the case of \code{noise = "quantile"} it is the mean of the eigenvalues smaller or equal to the \code{alpha}-quantile of the eigenvalues of the covariance matrix.

The augmentation estimator uses then the augmented components to measure the variation of the eigenvalues. For a more stable result it is recommened to repeat the augmentation process several times and Lue and Li (2021) recommend to use for \code{naug} approximately \code{p/5} or \code{p/10} where \code{p} is the number of columns of \code{X}. 

The augmented estimator for this purpose combines then the values of the 
scaled eigenvalues and the variation measured via augmentation. The main idea there is that for distinct eigenvales the variation of the eigenvectors
is small and for equal eigenvalues the corresponding eigenvectors have large variation.

The augmented estimate for k is the value where the measure takes its minimum and can be also visualized as a ladle.

For further details see Luo and Li (2021) and Radojicic et al. (2021).

}

\value{
A list of class ladle containing:
  \item{method}{the string PCA.}
  \item{k}{the estimated value of k.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
  \item{phin}{normalized eigenvalues of the covariance matrix.}
  \item{gn}{the main criterion for the augmented estimate - the sum of fn and phin. k is the value where gn takes its minimum}
  \item{lambda}{the eigenvalues of the covariance matrix.}
  \item{W}{the transformation matrix to the principal components.}
  \item{S}{data matrix with the centered principal components.}
  \item{MU}{the location of the data which was substracted before calculating the principal components.}
  \item{data.name}{the name of the data for which the augmented estimate was computed.}
  \item{sigma2}{the value used as noise variance when simulating the augmented components.}
}
\references{
\cite{Luo, W. and Li, B. (2021), On Order Determination by Predictor Augmentation, Biometrika, 108, 557--574. <doi:10.1093/biomet/asaa077>}

\cite{Radojicic, U., Lietzen, N., Nordhausen, K. and Virta, J. (2021), Dimension Estimation in Two-Dimensional PCA. In S. Loncaric, T. Petkovic and D. Petrinovic (editors) "Proceedings of the 12 International Symposium on Image and Signal Processing and Analysis (ISPA 2021)", 16--22. <doi:10.1109/ISPA52656.2021.9552114>}
}
\author{
Klaus Nordhausen
}
\seealso{
\code{\link{ladleplot}}, \code{\link{PCAladle}}
}
\examples{
n <- 1000
Y <- cbind(rnorm(n, sd=2), rnorm(n,sd=2), rnorm(n), rnorm(n), rnorm(n), rnorm(n))

testPCA <- PCAaug(Y) 
testPCA
summary(testPCA)
plot(testPCA)
ladleplot(testPCA)
ladleplot(testPCA, crit = "fn")
ladleplot(testPCA, crit = "lambda")
ladleplot(testPCA, crit = "phin")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
