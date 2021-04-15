\name{QPmin-package}
\alias{QPmin-package}
\docType{package}
\title{
\packageTitle{QPmin}
}
\description{
\packageDescription{QPmin}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{QPmin}
\packageIndices{QPmin}
Active Set Method (ASM) for minimizinig indefinite quadratic programs subject to
linear constraints: QPmin. The package also provide a set of additional functions
for the random generation of specific quadratic programs.
}
\author{
\packageAuthor{QPmin}

Maintainer: \packageMaintainer{QPmin}
}
\references{
``A new technique for generating quadratic programming test problems,'' Calamai P.H., L.N. Vicente, and J.J. Judice, Mathematical Programming 61 (1993), pp. 215-231.

``A General Quadratic Programming Algorithm,'' Fletcher, R., IMA Journal of Applied Mathematics 7 (1971), pp. 76-91.

``An Algorithm for Large-Scale Quadratic Programming,'' Gould, N.I.M., IMA Journal of Numerical Analysis 11 (1991), pp. 299-324. 

``A Stabilization of the Simplex Method,'' Bartels, R.H., Numerische Mathematik 16 (1971), pp. 414-434.
}

\keyword{ package }

\examples{
set.seed(2)
RP <- randomQP(8, "indefinite")
with(RP, QPmin(G, g, t(A), b, neq = 0, tol = 1e-8))
}
