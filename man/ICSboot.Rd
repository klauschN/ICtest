\name{ICSboot}
\alias{ICSboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boostrap-based Testing for the Number of Gaussian Components in NGCA Using Two Scatter Matrices
}
\description{
In independent components analysis (ICA) gaussian components are considered as uninteresting.
The function uses boostrappping tests, based on ICS using any combination of two scatter matrices, to decide if there are \code{p-k} gaussian components where \code{p} is the dimension of the data.
The function offers two different boostrapping strategies.
}
\usage{
ICSboot(X, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL, n.boot = 200, s.boot = "B1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{a numeric data matrix with p>1 columns.}
  \item{k}{the number of non-gaussian components under the null.}
  \item{S1}{name of the first scatter matrix function. Can only return a matrix. Default is \code{\link{cov}}}.
  \item{S2}{name of the second scatter matrix function. Can only return a matrix. Default is \code{\link[ICS]{cov4}}}
  \item{S1args}{list with optional additional arguments for \code{S1}.}
  \item{S2args}{list with optional additional arguments for \code{S2}.}
  \item{n.boot}{number of bootstrapping samples.}
  \item{s.boot}{bootstrapping strategy to be used. Possible values are \code{"B1"}, \code{"B2"}. See details for further information.}
}
\details{
While in \code{\link{FOBIasymp}} and \code{\link{FOBIboot}} the two scatters used are always \code{cov} and  \code{\link[ICS]{cov4}} this function can be used with any two scatter functions. In that case however the value of the Gaussian eigenvalues are in general not known and depend on the scatter functions used. Therefore the test uses as test statistic the \code{k} successive eigenvalues with the smallest variance. Which means the default here might differ from \code{\link{FOBIasymp}} and \code{\link{FOBIboot}}.

Given eigenvalues \eqn{d_1,...,d_p} the function thus orders the components in descending order according to the "variance" criterion .

Under the null it is then assumed that the first \code{k} interesting components are mutually independent and non-normal and the last \code{p-k} components are gaussian.

Let \eqn{d_1,...,d_p} be the ordered eigenvalues, \eqn{W} the correspondingly ordered unmixing matrix, \eqn{s_i = W (x_i-MU)}{s_i = W (X_i-MU)} the corresponding
source vectors which give the source matrix \eqn{S} which can be decomposed into \eqn{S_1} and \eqn{S_2} where \eqn{S_1} is the matrix with the \eqn{k} non-gaussian components
and \eqn{S_2} the matrix with the gaussian components (under the null).


Two possible bootstrap tests are provided for testing that the last \code{p-k} components are gaussian and independent from the first k components:

\enumerate{
\item \code{s.boot="B1"}: 
The first strategy has the followong steps:
\enumerate{
\item Take a bootstrap sample \eqn{S_1^*} of size \eqn{n} from \eqn{S_1}.
\item Take a bootstrap sample \eqn{S_2^*} consisting of a matrix with gaussian random variables having \eqn{cov(S_2)}.
\item Combine \eqn{S^*=(S_1^*, S_2^*)} and create \eqn{X^*= S^* W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}

Note that in this bootstrapping test the assumption of ''independent components'' is not used, it is only used that the last \eqn{p-k} components are gaussian and independent from the first \eqn{k} components. Therefore this strategy can be applied in an independent component analysis (ICA) framework
and in a non-gaussian components analysis (NGCA) framework.

\item \code{s.boot="B2"}: 
The second strategy has the following steps:
\enumerate{
\item Take a bootstrap sample \eqn{S_1^*} of size \eqn{n} from \eqn{S_1} where the subsampling is done separately for each independent component.
\item Take a bootstrap sample \eqn{S_2^*}  consisting of a matrix with gaussian random variables having \eqn{cov(S_2)}
\item Combine \eqn{S^*=(S_1^*, S_2^*)} and create \eqn{X^*= S^* W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}

This bootstrapping strategy assumes a full ICA model and cannot be used in an NGCA framework. Note that when the goal is to recover the non-gaussian independent components both scatters used must have the independence property.
}


}
\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the number of boostrapping samples used to obtain the p-value.}
  \item{method}{character string which test was performed and which scatters were used.}
  \item{data.name}{character string giving the name of the data.}
  \item{alternative}{character string specifying the alternative hypothesis.}
  \item{k}{the number or non-gaussian components used in the testing problem.}
  \item{W}{the transformation matrix to the independent components. Also known as unmixing matrix.}
  \item{S}{data matrix with the centered independent components.}
  \item{D}{the underlying eigenvalues.}
  \item{MU}{the location of the data which was substracted before calculating the independent components.}
  \item{s.boot}{character string which boostrapping strategy was used.}
}

\references{
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2022), Asymptotic and Bootstrap Tests for Subspace Dimension, Journal of Multivariate Analysis, 188, 104830. <doi:10.1016/j.jmva.2021.104830>.}

\cite{Nordhausen, K., Oja, H., Tyler, D.E. and Virta, J. (2017), Asymptotic and Bootstrap Tests for the Dimension of the Non-Gaussian Subspace, Signal Processing Letters, 24, 887--891. <doi:10.1109/LSP.2017.2696880>.}

\cite{Radojicic, U. and Nordhausen, K. (2020),  Non-Gaussian Component Analysis: Testing the Dimension of the Signal Subspace. In Maciak, M., Pestas, M. and Schindler, M. (editors) "Analytical Methods in Statistics. AMISTAT 2019", 101--123, Springer, Cham. <doi:10.1007/978-3-030-48814-7_6>.}
}

\author{
Klaus Nordhausen
}


\seealso{
\code{\link[ICS]{ics}}, \code{\link{FOBIboot}},  \code{\link{FOBIasymp}}
}

\examples{
n <- 750
S <- cbind(runif(n), rchisq(n, 2), rexp(n), rnorm(n), rnorm(n), rnorm(n))
A <- matrix(rnorm(36), ncol = 6)
X <- S \%*\% t(A)

# n.boot is small for demonstration purpose, should be larger
ICSboot(X, k=1, n.boot=20)

if(require("ICSNP")){

myTyl <- function(X,...) HR.Mest(X,...)$scatter
myT <- function(X,...) tM(X,...)$V

# n.boot is small for demonstration purpose, should be larger
ICSboot(X, k=3, S1=myT, S2=myTyl, s.boot = "B2", n.boot=20)
}
}

\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
