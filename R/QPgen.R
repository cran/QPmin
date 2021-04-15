QPgen.internal.scramble <- function(G, g, A, xstar, v, D) {
  n <- nrow(G)
  v <- v/sqrt(sum(v^2))
  H <- diag(n) - 2*v%*%t(v)
  M <- D%*%H
  Dinv <- diag(n)
  diag(Dinv) <- 1.0/diag(D)
  W <- H%*%Dinv
  Gout <- t(M)%*%G%*%M
  gout <- t(M)%*%g
  Aout <- A%*%M
  xstarout <- list()
  for(i in 1:length(xstar)) xstarout[[i]] <- W%*%xstar[[i]]
  return(list(G = Gout, g = gout, A = Aout, globals = xstarout))
}

QPgen.internal.concave <- function(m, thetas, alphas, betas, L) {  
  n <- 2*m
  gamma <- 3*m
  m0 <- sum(thetas == 0)
  m1 <- sum(thetas == 1)
  qs <- -(4^(1 - L))^(1 - thetas)
  ss <- (4^thetas)*qs
  xis <- (16^thetas)*qs
  G <- diag(c(qs,qs))
  g <- matrix(c(-ss, -ss), nrow = n, ncol = 1)
  xstar <- matrix(NA, nrow = n, ncol = 1)
  A <- matrix(0, nrow = gamma, ncol = n)
  b <- matrix(0, nrow = gamma, ncol = 1)
  for(l in 1:m) {
    A[3*l - 2, l] <- -alphas[l] 
    A[3*l - 2, m + l] <- -betas[l] 
    b[3*l - 2] <- -alphas[l] - betas[l] - alphas[l]*betas[l]
    A[3*l - 1, l] <- -1 
    A[3*l - 1, m + l] <- (betas[l] + 1) 
    b[3*l - 1] <- 0
    A[3*l, l] <- (alphas[l] + 1) 
    A[3*l, m + l] <- -1 
    b[3*l] <- 0
    if(thetas[l] == 0) {
      vals <- c(-0.5*alphas[l]*alphas[l]*4^(1 - L), 
                -0.5*betas[l]*betas[l]*4^(1 - L), 
                -4^(1 - L))
      iwin <- which.min(vals)
      if(iwin == 1) {
        xstar[l, 1] <- 1
        xstar[m + l, 1] <- 1 + alphas[l]
      } else if(iwin == 2) {
        xstar[l, 1] <- 1 + betas[l]
        xstar[m + l, 1] <- 1
      } else if(iwin == 3) {
        xstar[l, 1] <- xstar[m + l, 1] <- 0
      }
    } else {
      xstar[l, 1] <- xstar[m + l, 1] <- 0
    }
  }
  xstar <- list(xstar)
  optVal <- 0.5*t(xstar[[1]])%*%G%*%xstar[[1]] + t(g)%*%xstar[[1]]
  return(list(G = G, g = g, A = A, b = b, opt = optVal, globals = xstar))
}

QPgen.internal.bilinear <- function(m, alphas) {
  n <- 2*m #number of variables
  gamma <- 3*m #number of constraints
  g <- matrix(-1, nrow = n, ncol = 1)
  G <- matrix(0, n, n)
  diag(G[(m + 1):n, 1:m] ) <- 1
  diag(G[1:m, (m + 1):n]) <- 1
  A <- matrix(0, nrow = gamma, ncol = n)
  b <- matrix(NA, nrow = gamma, ncol = 1)
  mmin <- alphas < 0.5
  nmin <- sum(mmin)
  meq <- alphas == 0.5#there are 2^(mequal) global solutions
  neq <- sum(meq)
  mmaj <- alphas > 0.5
  nmaj <- sum(mmaj)
  alphas <- c(alphas[mmin], alphas[mmaj], alphas[meq]) #reordering makes it easier later on
  constant <- m
  for(l in 1:m) {
    A[3*l - 2, l] <- -alphas[l] 
    A[3*l - 2, m + l] <- -alphas[l] - 1 
    b[3*l - 2] <- -3*alphas[l] - 1
    A[3*l - 1, l] <- alphas[l] + 1 
    A[3*l - 1, m + l] <- alphas[l]
    b[3*l - 1] <- alphas[l] + 1
    A[3*l, l] <- -1
    A[3*l, m + l] <- 1
    b[3*l] <- -1
  }
  nGlo <- 2^neq
  xstar <- list()
  xstarbase <- matrix(NA, nrow = n, ncol = 1)
  if(nmin > 0) for(l in 1:nmin) {
    xstarbase[c(l, m + l), ] <- c(1.5, 0.5)
  }
  if(nmaj > 0) for(l in (nmin + 1):(nmin + nmaj)) {
    xstarbase[c(l, m + l), ] <- c(1 - alphas[l], 1 + alphas[l])
  }
  xstar[[1]] <- xstarbase
  addSolutions <- function(lev, theList) { 
    curLen <- length(theList)
    for(i in 1:curLen) {
      theList[[curLen + i]] <- theList[[i]]
      theList[[i]][c(lev, m + lev), ] <- c(1.5, 0.5)
      theList[[curLen + i]][c(lev, m + lev), ] <- c(0.5, 1.5)
    }
    if(lev < m) {
      theList = addSolutions(lev + 1, theList) 
      return(theList)
    }
    else return(theList)
  }
  if(neq > 0) xstar <- addSolutions(nmin + nmaj + 1, xstar)
  optVal <- 0.5*t(xstar[[1]])%*%G%*%xstar[[1]] + t(g)%*%xstar[[1]]
  return(list(G = G, g = g, A = A, b = b, opt = optVal, globals = xstar))
}

QPgen.internal.convex <- function(m, alphas, rhos, omegas) {
  n <- 2*m
  gamma <- 3*m
  thetas <- 1 - rhos*omegas
  xis <- 0.5*(1 + rhos)*(9^thetas)
  g <- matrix(c(-3^thetas, -rhos*3^thetas), nrow = n, ncol = 1)
  G <- diag(n)
  diag(G)[(m + 1):n] <- rhos
  A <- matrix(0, nrow = gamma, ncol = n)
  b <- matrix(NA, nrow = gamma, ncol = 1)
  xstar <- matrix(NA, nrow = n, ncol = 1)
  m01 <- thetas == 0 & rhos == 1
  m11 <- thetas == 1 & rhos == 1
  m10 <- thetas == 1 & rhos == 0
  for(l in 1:m) {
    A[3*l - 2, l] <- 3
    A[3*l - 2, m + l] <- 2
    b[3*l - 2] <- alphas[l]
    A[3*l - 1, l] <- 2
    A[3*l - 1, m + l] <- 3
    b[3*l - 1] <- alphas[l]
    A[3*l, l] <- -1
    A[3*l, m + l] <- -1
    b[3*l] <- -3
    if(m01[l]) {
      xstar[l, 1] <- xstar[m + l, 1] <- 0.2*alphas[l]
    } else if(m11[l]) {
      xstar[l, 1] <- xstar[m + l, 1] <- 1.5
    } else {
      xstar[l, 1] <- 9 - alphas[l]
      xstar[m + l, 1] <- alphas[l] - 6
    }
  }
  xstar <- list(xstar)
  optVal <- 0.5*t(xstar[[1]])%*%G%*%xstar[[1]] + t(g)%*%xstar[[1]]
  return(list(G = G, g = g, A = A, b = b, opt = optVal, globals = xstar))
}

randomQP <- function(n, type = c("convex", "concave", "indefinite")) {
  if(n%%2 != 0) {
    message("randomQP can only generate problems with an even number of variables. Setting n = n + 1.")
    n <- n + 1 
  }
  m <- n/2
  if(type == "convex") {
    alphas <- runif(m, min = 5, max = 7.4999)
    rhos <- round(runif(m))
    omegas <- round(runif(m))
    P <- QPgen.internal.convex(m, alphas, rhos, omegas)
  } else if(type == "indefinite") {
    numberOfHalfs <- ceiling(log(m))
    nmiss <- m - numberOfHalfs
    alphas <- c(runif(nmiss), rep(0.5, numberOfHalfs))
    P <- QPgen.internal.bilinear(m, alphas)  
  } else if (type == "concave") {
    thetas <- round(runif(m))
    draws <- runif(m)
    alphas <- betas <- c()
    for(l in 1:m) {
      if(draws[l] < 0.5) {
        alphas[l] <- 1.5
        betas[l] <- 2
      } else {
        alphas[l] <- 2
        betas[l] <- 1.5
      }
    }
    L <- ceiling(log(m))
    P <- QPgen.internal.concave(m, thetas, alphas, betas, L)
  }
  v <- matrix(0, nrow = 2*m, ncol = 1)
  draws <- runif(n)
  vals <- runif(n)
  v[draws > 0.35] <- vals[draws > 0.35] #put in some zeroes
  D <- diag(10*runif(n))
  PP <- QPgen.internal.scramble(P$G, P$g, P$A, P$globals, v, D)
  return(list(G = PP$G, g = PP$g, A = PP$A, b = P$b, opt = P$opt, solutions = PP$globals))
}









