\name{FOBIaug}
\alias{FOBIaug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Augmentation Estimate for the Number of Gaussian Components in ICA or NGCA using FOBI
}
\description{
The augmentation estimator uses the eigenvalues and eigenvectors of FOBI together with augmented Gaussian components to estimate the number of non-Gaussian components in ICA or NGCA.
}
\usage{
FOBIaug(X, naug = 1, nrep = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix.}
  \item{naug}{number of components to be augmented.}
  \item{nrep}{number of repetitions for the augmentation procedure.}
}
\details{
The model here assumes that in ICA or NGCA there are k non-gaussian components and p-k gaussian components.
The idea is then to decide which eigenvalues of the FOBI matrix differ from p+2. For this purpose the matrix
\eqn{M = (COV_4 - (p+2)I_p)^2}{M = (COV4 - (p+2)I)^2} is used, where \eqn{COV_4}{COV4} is the matrix of fourth moments of the
standardized data. Its eigenvalues are of the form \eqn{\lambda_1 \geq ... \geq \lambda_k > 0 =  ... = 0}{lambda_1 >= ... >= lambda_k > 0 =  ... = 0}.

The augmented estimator adds for that purpose \code{naug} standard Gaussian components to the data. As FOBI is affine equivariant the variance of the
augmented components does not matter and no estimation of a noise variance is needed. The augmented components belong then to the Gaussian part and the
augmentation estimator uses them to measure the variation of the eigenvectors. For a more stable result it is recommended to repeat the augmentation process
several times and Luo and Li (2021) recommend to use for \code{naug} approximately \code{p/5} or \code{p/10} where \code{p} is the number of columns of \code{X}.

The augmented estimator combines then the values of the scaled eigenvalues and the variation measured via augmentation. The main idea there is that for
distinct eigenvalues the variation of the eigenvectors is small and for equal eigenvalues the corresponding eigenvectors have large variation.

The augmented estimate for k is the value where the measure takes its minimum and can be also visualized as a ladle.

For further details see Luo and Li (2021).
}

\value{
A list of class ladle containing:
  \item{method}{the string FOBI.}
  \item{k}{the estimated number of non-gaussian components.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the augmented components for the different number of components.}
  \item{phin}{normalized eigenvalues of the FOBI matrix.}
  \item{gn}{the main criterion for the augmented estimate - the sum of fn and phin. k is the value where gn takes its minimum}
  \item{lambda}{the eigenvalues of the FOBI matrix.}
  \item{W}{the transformation matrix to the independent components. Also known as unmixing matrix.}
  \item{S}{data matrix with the centered independent components.}
  \item{MU}{the location of the data which was substracted before calculating the independent components.}
  \item{data.name}{the name of the data for which the augmented estimate was computed.}
}
\references{
\cite{Luo, W. and Li, B. (2021), On Order Determination by Predictor Augmentation, Biometrika, 108, 557--574. <doi:10.1093/biomet/asaa077>}
}
\author{
Klaus Nordhausen
}
\seealso{
\code{\link{ladleplot}}, \code{\link{FOBIladle}}, \code{\link{PCAaug}}, \code{\link{SIRaug}}
}
\examples{
n <- 1000
X <- cbind(rexp(n), rt(n,5), rnorm(n), rnorm(n), rnorm(n), rnorm(n))

test <- FOBIaug(X, nrep = 10)
test
summary(test)
plot(test)
ladleplot(test)
ladleplot(test, crit = "fn")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
