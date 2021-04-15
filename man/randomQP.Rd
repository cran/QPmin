\name{randomQP}
\alias{randomQP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
randomQP
}
\description{
Generates a non-separable random quadratic program of the specified type. 
}
\usage{
randomQP(n, type = c("convex", "concave", "indefinite"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{The random problem generated will have n variables and m = 3n/2 constraints. 
  Must be an even number.}
  \item{type}{Specifies the curvature of the objective function.}
}
\details{
The algorithm is based on Calamai, Vicente, and Judice (1993). It generates a 
random quadratic program with the following form
  \deqn{
  \min_x \frac{1}{2} x^T G x + x^T g 
  }
  \deqn{
   Ax \geq b
  }
}

\value{
  \item{G }{The quadratic component of the objective function. Must be symmetric.}
  \item{g}{The linear component of the objective function}
  \item{A}{The constraints coefficient matrix. This matrix has $n$ rows and $m$ columns.}
  \item{b}{The vector with the lower bounds on the constraints.}
  \item{opt}{An approximate expected value at the optimum solutions.}
  \item{solutions}{A list containing all of the global solutions to the problem.}
}
\references{
``A new technique for generating quadratic programming test problems,'' Calamai P.H., L.N. Vicente, and J.J. Judice, Mathematical Programming 61 (1993), pp. 215-231.
}
\author{
Andrea Giusto
}

\examples{
n <- 8

RP <- randomQP(n, "concave")

RP2 <- randomQP(n, "indefinite")

RP3 <- randomQP(n, "convex")

}
