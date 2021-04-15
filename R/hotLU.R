hotLU <- setClass(
  "hotLU",
  
  slots = c(
    tol = "numeric",
    O = "pMatrix",
    L = "dgCMatrix",
    U = "dgCMatrix",
    Q = "pMatrix",
    Lambdas = "list",
    Ps = "list"
  ),
  
#   prototype = list( L = NULL, U = NULL, Lambdas = list(), Ps = list()),
  
  validity = function(object) {
    nrow(object@L) -> rL
    ncol(object@L) -> cL
    nrow(object@U) -> rU
    ncol(object@U) -> cU
    if(rL != cL || rL != rU || rU != cU)
        return("L and U must be square and have matching dimensions.")
    return(TRUE)
  }
)

hotLU <- function(O = NULL, L = NULL, U = NULL, Q = NULL, tol = 1e-15) { 
  if(is.null(L)) L <- Matrix(0, nrow = 0, ncol = 0)
  if(is.null(U)) U <- Matrix(0, nrow = 0, ncol = 0)
  O <- as(O, "pMatrix")
  L <- drop0(as(L, "dgCMatrix"), tol)
  U <- drop0(as(U, "dgCMatrix"), tol)
  Q <- as(Q, "pMatrix")
  n <- nrow(L)
  if(is.null(O)) O <- as(1:n, "pMatrix")
  if(is.null(Q)) Q <- as(1:n, "pMatrix")
  new("hotLU", O = O, L = L, U = U, Q = Q, tol = tol)
}

setGeneric( name = "hotLUupdate", 
            def = function(or, p, spike) { 
              standardGeneric("hotLUupdate")
            }
          )
setMethod( f = "hotLUupdate",
           signature = "hotLU",
           definition = function(or, p, spike) {
              r <- length(or@Lambdas)
              n <- nrow(or@L)
              dc <- which(or@Q[p, ] == TRUE)
              Qupdate <- 1:n
              if(dc < n) {
                Qupdate[dc:(n-1)] <- Qupdate[(dc + 1):n]
                Qupdate[n] <- dc
              }
              Qupdate <- t(as(Qupdate, "pMatrix"))#pMatrix reasons row-wise, but my construction is column-major
              H <- or@U
              H[ ,dc] <- spike
              H <- H%*%Qupdate
              or@Q <- Qupdate%*%or@Q
              P <- as(1:n, "pMatrix")
              LambdaInv <- diag(n)
              if(dc < n) for(i in dc:(n - 1)) {
                lamInv <- as(diag(n), "dgCMatrix")
                curP <- as(1:n, "pMatrix")
                if(abs((P%*%H)[i, i]) < abs((P%*%H)[i + 1, i])) {
                  attr(P, "perm")[i] -> keep 
                  attr(P, "perm")[i] <- attr(P, "perm")[i + 1]
                  attr(P, "perm")[i + 1] <- keep
                  attr(curP, "perm")[i] -> keep
                  attr(curP, "perm")[i] <- attr(curP, "perm")[i + 1]
                  attr(curP, "perm")[i + 1] <- keep
                }
                factor <- (P%*%H)[i + 1, i]/(P%*%H)[i, i]
                lamInv[i + 1, i] <- factor
                LambdaInv <- as(LambdaInv%*%t(curP)%*%lamInv, "dgCMatrix")
                for(j in i:n)  
                  H[attr(P, "perm")[i + 1], j] <-  
                    H[attr(P, "perm")[i + 1], j] - factor*H[attr(P, "perm")[i], j] 
              }
              temp <- as(round(P%*%H, -floor(log10(.Machine$double.eps))), "dgCMatrix")
              or@U <- drop0(temp, or@tol)
              or@Lambdas[[r + 1]] <- drop0(P%*%LambdaInv, or@tol)
              or@Ps[[r + 1]] <- P
              return(or)
           }
         )

setGeneric( name = "hotLUsolve", 
            def = function(or, rhs) {
              standardGeneric("hotLUsolve")
            }
          )
setMethod( f = "hotLUsolve",
           signature = "hotLU",
           definition = function(or, rhs) {
              r <- length(or@Lambdas)
              spike <- forwardsolve(or@L, t(or@O)%*%rhs)
              if(r > 0) for(i in 1:r) 
                spike <- forwardsolve(or@Lambdas[[i]], or@Ps[[i]]%*%spike)
              sol <- or@Q%*%backsolve(or@U, spike)
              return(list(sol = sol, spike = spike))
           }
         )

setGeneric( name = "spike", 
            def = function(or, newCol) {
              standardGeneric("spike")
            }
          )
setMethod( f = "spike",
           signature = "hotLU",
           definition = function(or, newCol) {
              r <- length(or@Lambdas)
              spike <- forwardsolve(or@L, t(or@O)%*%newCol)
              if(r > 0) for(i in 1:r) 
                spike <- forwardsolve(or@Lambdas[[i]], or@Ps[[i]]%*%spike)
              return(spike)
           }
         )
## How to use
# nn <- 5
# set.seed(1)
# B <- Matrix(runif(nn*nn),nn,nn, sparse = TRUE) #Current basis
# LUD <- expand(lu(B))
# O <- t(LUD$P)
# L <- LUD$L
# U <- LUD$U
# Q <- t(LUD$Q)
# decompo <- hotLU(O, L, U, Q)
# 
# newCol1 <- runif(nn)
# spi <- spike(decompo, newCol1)
# spi - hotLUsolve(decompo, newCol1)$spike #good
# decompo <- hotLUupdate(decompo, p = 2, spike = spi)
# 
# newCol2 <- runif(nn)
# spi <- spike(decompo, newCol2)
# decompo <- hotLUupdate(decompo, p = 1, spike = spi)
# 
# newCol3 <- runif(nn)
# spi <- spike(decompo, newCol3)
# decompo <- hotLUupdate(decompo, p = 3, spike = spi)
# 
# B[ ,2] <- newCol1
# B[ ,1] <- newCol2
# B[ ,3] <- newCol3
# b <- 1:nn
# 
# solve(B, b)
# hotLUsolve(decompo, b)$sol #yay!


