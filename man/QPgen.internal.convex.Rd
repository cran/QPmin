\name{QPgen.internal.convex}
\alias{QPgen.internal.convex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QPgen.internal.convex
}
\description{
Generates a separable convex quadratic problem of the form
\deqn{
  \min_x \frac{1}{2} x^T G x + x^T g 
  }
  \deqn{
  Ax \geq b  
  }
}
\usage{
QPgen.internal.convex(m, alphas, rhos, omegas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Integer parameter controlling the number of variables (2m) and 
  constraints (3m) for the generated problem.}
  \item{alphas}{m parameters taking values between 5 and 7.5.}
  \item{rhos}{m parameters taking values in \{0, 1\}.}
  \item{omegas}{m parameters taking values in \{0, 1\}.}
}
\details{The problem has a unique global minimum and the constraints are 
linearly independent at all of the solutions.}
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
The function 'randomQP' uses 'QPgen.internal.convex' to construct non-separable 
convex problems. 
}


\examples{
m <- 3
alphas <-  runif(m, min = 5, max = 7.4999)
rhos <- round(runif(m))
omegas <- round(runif(m))
QPgen.internal.convex(m, alphas, rhos, omegas)
}
