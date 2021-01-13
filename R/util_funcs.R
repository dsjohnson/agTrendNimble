
iar.Q <- function (n, p){
  if (n < 2 * p) 
    stop("n must be >= 2*p\n")
  tmp1 <- out <- matrix(0, n, n)
  tmp2 <- ((-1)^c(0:p)) * choose(p, c(0:p))
  for (i in (p + 1):n) {
    tmp1[, i] <- c(rep(0, i - p - 1), tmp2, rep(0, n - i))
  }
  for (i in n:(p + 1)) {
    tmp4 <- tmp1[, c(1:n)[tmp1[i, ] != 0]]
    tmp5 <- tmp1[i, tmp1[i, ] != 0]
    tmp6 <- t(t(tmp4) * tmp5)
    tmp6[i, ] <- 0
    out[i, ] <- -rowSums(tmp6)
  }
  out[1:p, ] <- out[n:(n - p + 1), n:1]
  return(diag(apply(out, 1, sum)) - out)
}

iar_basis <- function(n, p, k){
  Q <- iar.Q(n,p)
  E <- eigen(Q)
  V <- E$vectors
  D <- 1/sqrt(zapsmall(E$values))
  D <- ifelse(!is.finite(D), 0, D)
  B <- V%*%diag(D)
  B <- B[,rev(order(D))]
  if(missing(k)) k<- n-p
  return(B[,1:k])
}


split_col <- function(M, idx){
  if(length(idx)!=ncol(M)) stop("length(idx) and ncol(M) are not equivalent!")
  u_idx <- 
  lapply(unique(idx), function(x, M,idx){M[,idx==x]}, M=M, idx=idx)
}




