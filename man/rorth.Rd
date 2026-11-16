\name{rorth}
\alias{rorth}
\title{
Random Orthogonal Matrix Uniformly Distributed with Respect to the Haar Measure
}
\description{
Creates a random orthogonal matrix which is uniformly distributed with respect to the Haar measure on the orthogonal group.
}
\usage{
rorth(k)
}
\arguments{
  \item{k}{the number of rows (and columns) of the orthogonal matrix.}
}

\details{
The function fills a \code{k} x \code{k} matrix \eqn{Z} with independent N(0,1) random variables and computes its QR decomposition \eqn{Z = QR} using \code{\link{qr}}. The orthogonal matrix \eqn{Q} is uniformly distributed with respect to the Haar measure only when the decomposition is made unique by requiring the diagonal elements of \eqn{R} to be positive. As \code{\link{qr}} does not impose this restriction, each column of \eqn{Q} is multiplied by the sign of the corresponding diagonal element of \eqn{R}.

Without this correction, \eqn{Q} is still orthogonal but not uniformly distributed. For example, the element \eqn{Q_{11}}{Q[1,1]} of \code{qr.Q(qr(Z))} is always negative, and the sign of its determinant depends only on \code{k}, whereas for a Haar distributed matrix both signs occur with probability 1/2. See the examples.

Such matrices are useful, for example, as random mixing matrices in simulation studies or for checking the affine equivariance of an estimator. In the package, \code{rorth} is used in the bootstrap test \code{\link{PCAboot}}.
}

\value{
An orthogonal \code{k} x \code{k} matrix.
}
\references{
\cite{Stewart, G.W. (1980), The efficient generation of random orthogonal matrices with an application to condition estimators, \emph{SIAM Journal on Numerical Analysis}, \bold{17}, 403--409. <doi:10.1137/0717034>.}

\cite{Mezzadri, F. (2007), How to generate random matrices from the classical compact groups, \emph{Notices of the American Mathematical Society}, \bold{54}, 592--604.}
}
\author{
Klaus Nordhausen
}

\examples{
Orth <- rorth(4)
crossprod(Orth)
tcrossprod(Orth)

# comparison with the uncorrected QR decomposition: for a Haar
# distributed matrix both proportions should be close to 1/2
naiveOrth <- function(k) qr.Q(qr(matrix(rnorm(k * k), k)))
sapply(list(qr = naiveOrth, rorth = rorth), function(f) {
  Q <- replicate(500, f(4), simplify = FALSE)
  c("P(Q[1,1] < 0)" = mean(sapply(Q, function(q) q[1, 1] < 0)),
    "P(det < 0)" = mean(sapply(Q, det) < 0))
})
}
\keyword{multivariate}
\keyword{datagen}
