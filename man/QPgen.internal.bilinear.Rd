\name{QPgen.internal.bilinear}
\alias{QPgen.internal.bilinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QPgen.internal.bilinear
}
\description{
Generates a separable bilinear problem of the form
\deqn{
  \min_x \frac{1}{2} x^T G x + x^T g 
  }
  \deqn{
  Ax \geq b  
  }
}
\usage{
QPgen.internal.bilinear(m, alphas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Integer parameter controlling the number of variables (2m) and 
  constraints (3m) for the generated problem.}
  \item{alphas}{m positive parameters.}
}
\details{The problem is an indefinite problem with 2^m local minima of which 2^n 
are global. Here, n is equal to the number of alphas exactly equal to 0.5. The
 constraints are guaranteed independent only at each solution, but not generally
everywhere in the feasible region.}
\value{
\item{G}{The quadratic component of the objective function.}
  \item{g}{The linear component of the objective function}
  \item{A}{The constraints coefficient matrix. This matrix has 3m rows and 2m columns.}
  \item{b}{The vector with the lower bounds on the constraints.}
  \item{opt}{An approximate expected value at the optimum solutions.}
  \item{globals}{A list containing all of the global solutions to the problem.}
}
\references{
``A new technique for generating quadratic programming test problems,'' Calamai P.H., L.N. Vicente, and J.J. Judice, Mathematical Programming 61 (1993), pp. 215-231.
}
\author{
Andrea Giusto
}
\note{
The function 'randomQP' uses 'QPgen.internal.bilinear' to construct non-separable 
indefinite problems. The technique used to conceal the separability of the problem
also eliminates bi-linearity.
}


\examples{
m <- 3
nhalves <- ceiling(log(m))
nmiss <- m - nhalves
alphas <- c(runif(nmiss), rep(0.5, nhalves))
QPgen.internal.bilinear(m, alphas)
}
