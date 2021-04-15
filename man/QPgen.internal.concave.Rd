\name{QPgen.internal.concave}
\alias{QPgen.internal.concave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QPgen.internal.concave
}
\description{
Generates a separable strictly concave quadratic problem of the form
\deqn{
  \min_x \frac{1}{2} x^T G x + x^T g 
  }
  \deqn{
  Ax \geq b  
  }
}
\usage{
QPgen.internal.concave(m, thetas, alphas, betas, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Integer parameter controlling the number of variables (2m) and 
  constraints (3m) for the generated problem.}
  \item{thetas}{m binary \{0, 1\} parameters.}
  \item{alphas}{m parameters taking on either of the two values \{1.5, 2\}}
  \item{betas}{m parameters taking on either of the two value \{1.5, 2\}. Each entry
  of betas must be different from each entry of alphas.}
  \item{L}{A positive parameter}
}
\details{Denoting m0 the number of thetas == 0, the generated problem has a unique
global minimum and 3^m0 local ones. All of the local minima have distinct function 
values. While the constraints are guaranteed to be indipendent at each solution, they 
are not guaranteed independent everywhere, for this reason too large problems tend to
not be solvable with QPmin, which relies on a basic version of the simplex method for
linear programming to find an initial feasible point. A large enough problem will 
generally halt the execution due to linear dependence in the basis constraints.}
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
The function 'randomQP' uses 'QPgen.internal.concave' to construct non-separable 
concave problems. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
m <- 3
thetas <- round(runif(m))
draws <- runif(m)
belowhalf <- draws < 0.5 
alphas <- betas <- c()
alphas[belowhalf] <- 1.5
alphas[!belowhalf] <- 2
betas[belowhalf] <- 2
betas[!belowhalf] <- 1.5
L <- ceiling(log(m))
QPgen.internal.concave(m, thetas, alphas, betas, L)
}
