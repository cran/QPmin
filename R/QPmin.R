QPmin <- function(G, g, A, b, neq, 
                        Lb = NULL, Ub = NULL, tol = 1e-6,
                        initialPoint = NULL, existingLU = NULL,
                        returnLU = FALSE) {
  G <- as(G, "dgCMatrix")
  n <- nrow(G)
  g <- as(g, "sparseVector")
  borig <- b <- as(b, "sparseVector")
  A <- as(A, "dgCMatrix")
  Atriplet <- as(A, "dgTMatrix")
  Atriplet <- t(Atriplet)
  m <- ncol(A)
  rmost <- ceiling(1.1*n)
  if(is.null(Lb)) Lb <- rep(-Inf, n)
  if(is.null(Ub)) Ub <- rep(Inf, n)
  temp <- countBounds(n, Lb, Ub)
  nlact <- temp$nlact; lact <- temp$lact
  nuact <- temp$nuact; uact <- temp$uact
  b <- c(b, temp$rhs)
  nctot <- m + nlact + nuact
  if(!is.null(initialPoint)) {
    geq0 <- A%*%initialPoint - borig
    if(any(geq0 < -tol)) {
      message("The initial point provided is not feasible.")
      return("unfeasible")
    }
    activeLBounds <- which(abs(initialPoint - Lb) < tol)
    activeUBounds <- which(abs(initialPoint - Ub) < tol)
    feas <- list(x = initialPoint, binding = which(abs(geq0) < tol), 
                 unfeasibles = NULL, activeLBounds = activeLBounds, 
                 activeUBounds = activeUBounds)
  } else if(n - nctot > 0) {  # the problem does not have enough constraints to use the 
    Laug <- findArbitraryBounds(n, m, A, lact, Lb, Ub, b, tol) 
    feas <- feasibility(Atriplet, borig, neq, Laug, Ub, tol)
  } else feas <- feasibility(Atriplet, borig, neq, Lb, Ub, tol)
  if(length(feas) == 1) {
    if(feas == "unbounded") return("unbounded") else return("singular")
  } else if(length(feas$unfeasibles) == 0) {
    ##Resizing an existing sparse matrix is not straightforward; resize and fill 
    ##with additional 1's and -1's for lower/upper constraints
    attr(A, "Dim") <- as.integer(c(n, nctot))
    lastn <- attr(A, "p")
    lastn <- as.integer(lastn[length(lastn)])
    if(nlact > 0) for(j in 1:nlact) {
      attr(A, "p") <- c(attr(A, "p"), as.integer(lastn + j))
      attr(A, "i") <- c(attr(A, "i"), as.integer(which(lact)[j] - 1))
      attr(A, "x") <- c(attr(A, "x"), 1)
    }
    if(nuact > 0) for(j in 1:nuact) {
      attr(A, "p") <- c(attr(A, "p"), as.integer(lastn + nlact + j))
      attr(A, "i") <- c(attr(A, "i"), as.integer(which(uact)[j] - 1))
      attr(A, "x") <- c(attr(A, "x"), -1)
    }
    x <-as(feas$x, "sparseVector")
    active <- feas$binding
    nn <- length(feas$activeLBounds)
    if(nlact > 0 && nn > 0)
      for(i in 1:nn) 
        if(lact[feas$activeLBounds[i]]) {
          idx <- 0
          for(ii in 1:feas$activeLBounds[i]) {
            if(lact[ii]) idx <- idx + 1
          }
          active <- c(active, m + idx)
        }
    nn <- length(feas$activeUBounds)
    if(nuact > 0 && nn > 0) 
      for(i in 1:nn) 
        if(uact[feas$activeUBounds[i]]) {
          idx <- 0
          for(ii in 1:feas$activeUBounds[i]) {
            if(lact[ii]) idx <- idx + 1
          }
          active <- c(active, m + nlact + idx)
        }
  } else {return("unfeasible")}
  nactive <- length(active)
  if(nactive > n) { 
    active <- active[1:n]
    nactive <- n
  }
  inactive <- setdiff(1:nctot, active)
  rhs <- as(-G%*%x - g, "sparseVector")
  attr(rhs, "length") <- as.integer(n + nctot)
  doit <- FALSE
  if(is.null(existingLU)) {
    doit <- TRUE
  } else {
    if(length(existingLU@Ps) > rmost) doit <- TRUE
  }
  if(doit) {
    Btop <- cbind(G, A[,active], matrix(0, nrow = n, ncol = nctot - nactive))
    Bmid <- cbind(t(A[,active]), matrix(0, nrow = nactive, ncol = nctot))
    Bbot <- cbind(t(A[,inactive]), matrix(0, nrow = nctot - nactive, ncol = nactive), diag(nctot - nactive))
    B <- rbind(Btop, Bmid, Bbot)
    LUD <- expand(lu(B))
    O <- t(LUD$P)
    L <- LUD$L
    U <- LUD$U
    Q <- t(LUD$Q)
    existingLU <- hotLU(O, L, U, Q, tol)
  }
  lambdas <- hotLUsolve(existingLU, rhs)$sol[(n + 1):(n + nactive)]
  if(nctot > neq) {
    mask <- (neq + 1):nactive
    maxL <- max(lambdas[mask]) 
  } else {
    maxL <- max(lambdas)
  }
  while(maxL > -tol) {
    jj <- which(lambdas[mask] == maxL)[1] + neq
    P <- 1:(n + nctot)
    P[n + jj] -> keep
    P[n + jj] <- P[n + nactive]
    P[n + nactive] <- keep
    P <- as(P, "pMatrix")
    existingLU@O <- existingLU@O%*%P
    existingLU@Q <- existingLU@Q%*%P
    lambdas[jj] -> keep
    lambdas[jj] <- lambdas[nactive]
    lambdas[nactive] <- keep
    active[jj] -> keep
    active[jj] <- active[nactive]
    active[nactive] <- keep
    rhs <- sparseVector(x = 1, i = n + nactive, length = n + nctot)
    theSolution <- hotLUsolve(existingLU, rhs)
    foo <- theSolution$sol
    p <- as(foo[1:n], "sparseVector")
    if(nactive > 1) q <- as(foo[(n + 1):(n + nactive - 1)], "sparseVector")
    mu <- foo[n + nactive]
    halts <- t(p)%*%A
    encountered <- which(halts <= - tol)
    alphas <- (b[encountered] - t(x)%*%A[,encountered])/halts[encountered]
    foundBlockers <- length(alphas) > 0
    if(abs(mu) < tol) alphaHat <- Inf else alphaHat <- as.numeric(-lambdas[nactive]/mu) 
    alphaBar <- 0
    if(foundBlockers) {
      alphaBar <- min(alphas)
      ll <- which(alphas == alphaBar)[1]
      l <- encountered[ll]
    }
    if((mu < -tol && !foundBlockers) || 
       (mu < -tol && alphaHat <= alphaBar)) {
#       CASE 1: drop active[nactive]
       j <- active[nactive]
       inactive <- c(j, inactive)
       active <- active[-nactive]
       lambdas <- lambdas[-nactive]
       nactive <- nactive - 1
       x <- x + alphaHat*p
       if(nactive > 0) lambdas <- lambdas + alphaHat*q
       spike <- theSolution$spike
       coltoupdate <- c(n + nactive + 1)
    } else if(mu > -tol && !foundBlockers) {
#       CASE 2
        return("unbounded")
    } else if(abs(mu) < tol && abs(maxL) < tol) {
#       CASE 3: a saddle has been encountered...TODO: finish this case
      if(any(alphaBar*q > tol)) {
        return("finish the get-out-of-a-saddle case")
      } else return(list(x = x, binding = active, lambdas = lambdas, ## all points 
                    saddle = x + alphaBar*p)) ## between x and x + alphaBar*p are solutions     
    } else {
      ll <- which(inactive == l)
      rhs <- as(A[,l], "sparseVector")
      attr(rhs, "length") <- n + nctot
      theNewSolution <- hotLUsolve(existingLU, rhs)
      foo <- theNewSolution$sol
      u <- as(foo[1:n], "sparseVector")
      if (nactive > 1) v <- as(foo[(n + 1) : (n + nactive - 1)], "sparseVector")
      theta <- foo[n + nactive]
      if(all(abs(u) < tol)) {
#         CASE 4: drop active[nactive] and add l = inactive[ll]
        P <- 1:(n + nctot)
        P[n + nactive] <- n + nactive + ll
        P[n + nactive + ll] <- n + nactive
        P <- as(P, "pMatrix")
        existingLU@O <- existingLU@O%*%P
        existingLU@Q <- existingLU@Q%*%P
        spike <- theNewSolution$spike
        coltoupdate <- c(n + nactive, n + nactive + ll)
        inactive[ll] <- active[nactive]
        active[nactive] <- l
        b[inactive[ll]] <- b[inactive[ll]] + runif(1, -tol, tol)
        x <- x + alphaBar*p        
        lambdas[nactive] <- (lambdas[nactive] + alphaBar*mu)/theta
        if(nactive > 1) lambdas[-nactive] <- lambdas[-nactive] + 
                                              alphaBar*q - lambdas[nactive]*v
      } else {
#         CASE 5: add l as active constraint
        P <- 1:(n + nctot)
        P[n + nactive + 1] -> keep
        P[n + nactive + 1] <- P[n + nactive + ll]
        P[n + nactive + ll] <- keep
        P <- as(P, "pMatrix")
        existingLU@O <- existingLU@O%*%P
        existingLU@Q <- existingLU@Q%*%P
        spike <- theNewSolution$spike
        coltoupdate <- n + nactive + 1
        lambdas[nactive] <- lambdas[nactive] + alphaBar*mu
        if(nactive > 1) lambdas[-nactive] <- lambdas[-nactive] + alphaBar*q
        active <- c(active, l) 
        inactive[ll] <- inactive[1]
        inactive <- inactive[-1]
        nactive <- nactive + 1
        lambdas <- c(lambdas, 0)
        x <- x + alphaBar*p
      }
    } 
    maxL <- max(lambdas)
    mask <- (neq + 1):nactive
    if(length(existingLU@Ps) <= rmost) {
      existingLU <- hotLUupdate(existingLU, coltoupdate[1], spike)
      if(length(coltoupdate) > 1) {
        rhs <- sparseVector(x = 1.0, i = coltoupdate[2], length = n + nctot)
        spiek <- spike(existingLU, rhs)
        existingLU <- hotLUupdate(existingLU, coltoupdate[2], spiek)
      }
    } else { # redo the LU factorization from scratch
      Btop <- cbind(G, A[,active], matrix(0, nrow = n, ncol = nctot - nactive))
      Bmid <- cbind(t(A[,active]), matrix(0, nrow = nactive, ncol = nctot))
      Bbot <- cbind(t(A[,inactive]), matrix(0, nrow = nctot - nactive, ncol = nactive), diag(nctot - nactive))
      B <- rbind(Btop, Bmid, Bbot)
      LUD <-  expand(lu(B))
      O <- t(LUD$P)
      L <- LUD$L
      U <- LUD$U
      Q <- t(LUD$Q)
      existingLU <- hotLU(O, L, U, Q, tol)
    }
  }
  res <- list(x = x, binding = active, lambdas = lambdas)
  if(returnLU) res[["warmStart"]] <- existingLU
  return(res)
}
