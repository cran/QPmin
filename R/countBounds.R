countBounds <- function(n, Lb, Ub) {
  nlact <- nuact <- 0
  lact <- uact <- c()
  rhs <- c()
  for(i in 1:n) {
    if(Lb[i] > -Inf) { 
      nlact <- nlact + 1
      lact[i] <- TRUE
      rhs <- c(rhs, Lb[i]);
    } else lact[i] <- FALSE
  } 
  for(i in 1:n) {
    if(Ub[i] < Inf) {
      nuact <- nuact + 1
      uact[i] <- TRUE
      rhs <- c(rhs, -Ub[i])
    } else uact[i] <- FALSE
  }
  if(is.null(rhs)) rhs <- NA
  return(list(lact = lact, nlact = nlact, uact = uact, nuact = nuact, 
              rhs = matrix(rhs, nlact + nuact, 1)))
}

findArbitraryBounds <- function(n, m, A, lact, Lb, Ub, b, tol) {
  hiddenLbounds <- rep(FALSE, n) # simplex-method to find the starting point.
  hiddenU <- rep(Inf, n)
  for(j in 1:m) { # Imposing additional lower bounds on all variables that 
    nnz <- 0      # are unbounded from below, and ensuring that these are 
    i <- 1        # compatible with hidden and existing upper bounds.
    inz <- 0
    while(nnz <= 1 && i <= n) { 
      if(abs(A[i, j]) > tol) {
          inz <- i
          nnz <- nnz + 1
        }
        i <- i + 1
    }
    if(nnz <= 1) {
      if(A[inz, j] > 0) { 
        hiddenLbounds[inz] <- TRUE
      } else {
        hiddenU[inz] <- -b[j]
      }
    }
  }
  lactAug <- lact | hiddenLbounds
  Uaug <- apply(cbind(Ub, hiddenU), 1, min)
  Laug <- Lb
  i <- 1
  while(i <= n) {
    if(!lactAug[i]) {
      if(0 < Uaug[i]) {
        Laug[i] <- 0
      } else if(-1 < Uaug[i]) { #in case upper bound is zero
        Laug[i] <- -1 
      } else Laug[i] <- Uaug[i] - 1
    }
    i <- i + 1
  }
  return(Laug)
}



