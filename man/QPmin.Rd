\name{QPmin}
\alias{QPmin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QPmin
}
\description{
Active set method for the solution of the quadratic program
\deqn{
  \min_x \frac{1}{2} x^T G x + x^T g 
  }
  \deqn{
  A^Tx \geq b  
  }
  \deqn{
  Lb \leq x \leq Ub
  }
  where the first neq rows of the system \eqn{A^Tx \geq b} are strict equalities.
}
\usage{
QPmin(G, g, A, b, neq, Lb = NULL, Ub = NULL, tol = 1e-06, initialPoint = NULL, 
      existingLU = NULL, returnLU = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{The quadratic term of the objective function. Must be symmetric.}
  \item{g}{The linear term of the objective function.}
  \item{A}{The constraints coefficient matrix. This matrix has as many rows as the 
  dimensions of G and g and an arbitrary number of columns. If the problem has 
  no linear constraints (but it still may have simple bound constraints), it is 
  still necessary to pass a matrix with the correct number of rows and zero columns.}
  \item{b}{The lower bounds on the constraints.}
  \item{neq}{The number of rows of A to be interpreted as strict equalities.}
  \item{Lb}{An array with the lower bounds on the variables. If there are no 
  lower bounds, this matrix may be omitted, but if otherwise specified, it must 
  contain an entry for every variable of the problem in the corresponding 
  position. The value -Inf can be used to specify that a variable is 
  unconstrained from below.}
  \item{Ub}{An array with the upper bounds on the variables. The rules about Lb 
  apply to this parameter, too.}
  \item{tol}{A tolerance parameter for the evaluation of the convergence criterion.}
  \item{initialPoint}{Intended as a "hot start" facility, when the user has a feasible
  starting point for which the second order conditions are satisfied. The code checks 
  if initialPoint is feasible and it terminates if it isn't. However, there are no 
  checks to insure thta the second order conditions are satisfied. See Fletcher (1971).}
  \item{existingLU}{Intended as a "warm start" facility, it allows the user to avoid the 
  re-factorization of a valid basis matrix into its L and U components. See Fletcher 
  (1971) and Bartels (1969).}
  \item{returnLU}{Logical value controlling if the updated LU decomposition of the 
  basis matrix should be returned for later re-use as part of a warm-start strategy.}
}

\value{
 \item{x}{The (possibly local) solution found by the algorithm. }
 \item{active}{A list with the active constraints. The first neq constraints are always
  binding. The numbers above m relate to the lower and upper bounds respectively as 
  specified by the user. For example in a case with 8 variables 5 constraints 
  (2 of which are equalities), with lower bounds on the third and fifth 
  variables, and an upper bound on the first variable,  an active set of 
  \{1, 2, 4, 6, 8 \} indicates that the first, second and fourth constraints are 
  binding, as well as both the lower bounds on variable 3 and upper bound on 
  variable 1.}
  \item{lambdas}{The Lagrange multipliers associated with each one of the binding 
  constraints}
  \item{warmStart}{The existing LU decoposition in its product form.}
}
\references{
``A General Quadratic Programming Algorithm,'' Fletcher, R., IMA Journal of Applied Mathematics 7 (1971), pp. 76-91.

``An Algorithm for Large-Scale Quadratic Programming,'' Gould, N.I.M., IMA Journal of Numerical Analysis 11 (1991), pp. 299-324. 

``A Stabilization of the Simplex Method,'' Bartels, R.H., Numerische Mathematik 16 (1971), pp. 414-434.
}
\author{
Andrea Giusto
}

\examples{
G <- diag(5)
g <- c(0,-6,-6,-12,-9)
A <- matrix(c(2, rep(0,3), -1, 5, 0, -3, 0, -1, 0, -1, 0, -3, 0), 5, 3)
b <- c(0, 0, 0)
Lb <- c(-Inf, -Inf, 0, 0, 0)
Ub <- c(4, 8, Inf, Inf, Inf)
neq <- 0
QPmin(G, g, A, b, neq = neq, Lb, Ub, tol = 1e-6)

set.seed(2)
RP <- randomQP(8, "indefinite")
Lb <- rep(-Inf, 8)
Ub <- -Lb 
with(RP, QPmin(G, g, t(A), b, neq = 0, Lb, Ub, tol = 1e-06))


}
