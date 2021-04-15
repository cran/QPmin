# Simple implementation of the simplex algorithm for the identification of a 
# feasible x given a set of constraints in the form Ax >= b.  The matrix A is 
# assumed row-major and the first neq rows are strict equalities. On return  
# the most feasible x is returned in with the list of binding constraints and 
# active bounds. The bounds are reported by variable (e.g. 2 in the 
# activeUBounds list indicates that variable 2 is at the lower bound). If the 
# list of infeasibles is not empty, the problem is not feasible.
# At is assumed passed as a "dgTMatrix" which makes it easy to construct the 
# constraint matrix including the slack and artificial variables.
feasibility <- function(Atriplet, b, neq, Lb, Ub, tol = 1e-6) {
  n <- ncol(Atriplet)
  m <- nrow(Atriplet)
  lact <- Lb > -Inf
  uact <- Ub < Inf
  nlact <- sum(lact)
  nuact <- sum(uact)
  nctot <- m + nlact + nuact
  xm <- xp <- matrix(0, n, 1)
  ##constructing an empty sparse matrix 
  ri <- ci <- vals <- c()
#   AA <- sparseMatrix(0, nctot, 2*n + 2*nctot - neq)
  ##Original coefficients first: AA[1:m, 1:(2*n)] <- cbind(Atriplet, -Atriplet) 
  nnz <- 0
  if(m > 0) {
    ri <- c(attr(Atriplet, "i"), attr(Atriplet, "i")) + 1 
    ci <- c(attr(Atriplet, "j"), attr(Atriplet, "j") + n) + 1
    vals <- c(attr(Atriplet, "x"), -attr(Atriplet, "x"))
    nnz <- length(vals)
  }
  nnz <- nnz + 1
  ##Slack variables for inequalities: for(j in 1:(m - neq)) AA[neq + j, j + 2*n] <- -1
  if(m - neq > 0) for(j in 1:(m - neq)) {
    ri[nnz] <- neq + j
    ci[nnz] <- 2*n + j
    vals[nnz] <- -1
    nnz <- nnz + 1
  }
  ##Artificial for original constraints: for(j in 1:m) AA[j, j + 2*n + nctot - neq] <- sign(b[j] + tol)
  if(m > 0) for(j in 1:m) {
    ri[nnz] <- j
    ci[nnz] <- 2*n + nctot - neq + j
    vals[nnz] <- sign(b[j] + tol)
    nnz <- nnz + 1
  }
  ##Lower bounds if present AA[(m + 1):(m + nlact), 1:(2*n)] <- cbind(Ibar, -Ibar)
  ##                        for(j in 1:nlact) AA[j + m, j + 2*n + m - neq] <- -1
  ##                        for(j in 1:nlact) AA[j + m, j + 2*n + nctot - neq + m] <- sign(b[j + m] + tol)
  if(nlact > 0) for(j in 1:nlact) {
    ri[nnz] <- m + j
    whichVar <- which(lact)[j] 
    ci[nnz] <-  whichVar
    vals[nnz] <- 1
    nnz <- nnz + 1
    ri[nnz] <- m + j
    ci[nnz] <- n + whichVar 
    vals[nnz] <- -1
    nnz <- nnz + 1
    ##Add the rhs coefficient
    b <- c(b, Lb[whichVar])
    ##Slack and Artificial variables for existing lower bounds
    ri[nnz] <- m + j
    ci[nnz] <- 2*n + m - neq + j
    vals[nnz] <- -1
    nnz <- nnz + 1
    ri[nnz] <- m + j
    ci[nnz] <- 2*n + nctot - neq + m + j
    vals[nnz] <- sign(Lb[whichVar] + tol)
    nnz <- nnz + 1
  }
  ##Upper bounds if present AA[(m + nlact + 1):(nctot), 1:(2*n)] <- cbind(-barI, barI)
  ## for(j in 1:nuact) AA[j + m + nlact, j + 2*n + m - neq + nlact] <- -1
  ##                   AA[j + m + nlact, j + 2*n + nctot - neq + nlact + m] <- sign(b[j + m + nlact] + tol)
  if(nuact > 0) {
    for(j in 1:nuact) {
      ri[nnz] <- m + nlact + j
      whichVar <- which(uact)[j]  
      ci[nnz] <- whichVar 
      vals[nnz] <- -1
      nnz <- nnz + 1
      ri[nnz] <- m + nlact + j
      ci[nnz] <- n + whichVar
      vals[nnz] <- 1
      nnz <- nnz + 1
      ##Add the rhs coefficient
      b <- c(b, -Ub[whichVar])
      ####Slack and Artificial variables for existing upper bounds
      ri[nnz] <- m + nlact + j
      ci[nnz] <- 2*n + m - neq + nlact + j
      vals[nnz] <- - 1
      nnz <- nnz + 1
      ri[nnz] <- m + nlact + j
      ci[nnz] <- 2*(n + m + nlact) + nuact - neq + j
      vals[nnz] <- sign(-Ub[whichVar] + tol)
      nnz <- nnz + 1
    }
  }
  AA <- sparseMatrix(i = ri, j = ci, x = vals)
  ibeg <- 2*n + nctot - neq + 1
  c <- sparseVector(x = rep(1, nctot), i = ibeg:ncol(AA), 
                    length = ncol(AA))
  Ba <- ibeg:ncol(AA)
  Na <- setdiff(1:ncol(AA), Ba)
  Binv <- AA[,ibeg:ncol(AA)]
  stop <- FALSE
  N <- AA[,Na,drop=FALSE]
  while(!stop) {
    xb <- as(Binv%*%b, "sparseVector") ##NOTE: sparseVector is base1 apparently :/
    lambda <- as(c[Ba]%*%Binv, "sparseVector")
    sN <- as(c[Na] - lambda%*%N, "sparseVector")
    if(all(sN >= -tol)) 
    {
      stop <- TRUE
    } else {
      qq <- which.min(attr(sN, "x"))
      qq <- attr(sN, "i")[qq]
      q <- Na[qq]
      wtil <- as(Binv%*%AA[,q], "sparseVector")
      if(all(wtil < tol)) {
        return("unbounded")
      }
      d <- wtil
      d[d < 0] <- 0
      temp <- xb/d ##NOTE: apparently sparseVector/sparseVector -> numeric :/
      pps <- which(temp  == min(temp, na.rm = TRUE))
      maxwtil <- max(wtil[pps])
      pp <- pps[which(wtil[pps] == maxwtil)]
      if(abs(wtil[pp]) < tol) return("singular")
      p <- Ba[pp]
      Ba[pp] <- q
      Na[qq] <- p
      factor <- as.numeric(wtil[pp])
      wtil[pp] <- wtil[pp] - 1
      wtil <- -wtil/factor
      ri <- ci <- vals <- c()
      if(pp > 1) {
        ci <- ri <- 1:(pp-1)
        vals <- rep(0, (pp - 1))
      }
      ri <- c(ri, attr(wtil, "i"))
      ci <- c(ci, rep(pp, length(attr(wtil, "x"))))
      vals <- c(vals, attr(wtil, "x"))
      if(pp < nctot) {
        ri <- c(ri, (pp + 1):nctot)
        ci <- c(ci, (pp + 1):nctot)
        vals <- c(vals, rep(0, nctot - pp))
      }
      Etil <- sparseMatrix(i = ri, j = ci, x = vals, dims = c(nctot, nctot))
      diag(Etil) <- diag(Etil) + 1
      Binv <- Etil%*%Binv
      N[,qq] <- AA[,p]
#       print(Ba)
#       print(round(Binv, 4))
#       print(stop)
    }
  }
  sol <- rep(0, ncol(AA))
  for(i in 1:nctot) sol[Ba[i]] = xb[i]
  pl <- 1:n
  mn <- (n+1):(2*n)
  x <- sol[pl] - sol[mn]
  if(neq > 0) binders <- 1:neq else binders <- c()
  if(m - neq > 0) binders <- c(binders, neq + which(sol[(2*n + 1):(2*n + m - neq)] < tol))
  activeUBounds <- activeLBounds <- NULL
  if(nlact > 0) {
    activeLBounds <- which(sol[(2*n + m - neq + 1):(2*n + m - neq + nlact)] < tol)
    activeLBounds <- which(lact == TRUE)[activeLBounds]
  }
  if(nuact > 0) {
    activeUBounds <- which(sol[(2*n + m - neq + nlact + 1):(2*n + nctot - neq)] < tol)
    activeUBounds <- which(uact == TRUE)[activeUBounds]
  }
  unfeas <- which(sol[(2*n + nctot - neq + 1):length(sol)] > tol)
  return(list(x = x, binding = binders, activeLBounds = activeLBounds, 
              activeUBounds =  activeUBounds, unfeasibles = unfeas) )
}

