#########################################
#### SLP algorithm in Overton (1992) ####
#########################################
# pp. 104-105 of Overton (1992) with full linear programming instead of PLP algorithm
# TO-DO: use efficient eigen-decomposition due to dpr1 structure

#############################
#### svec transformation ####
#############################
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

svec <- function(X){
  # symmetric X of size n times n
  ind <- as.vector(row(X)> col(X))
  ind2 <- as.vector(row(X)>= col(X))
  X[ind] <- X[ind] * sqrt(2)
  return(X[ind2])
}

svec.inv <- function(x){
  nn <- length(x)
  n <- (-1 + sqrt(1+4 *2*nn))/2
  X <- matrix(nr=n, nc=n)
  ind <- as.vector(row(X)> col(X))
  ind2 <- as.vector(row(X)>= col(X))
  ind3 <- as.vector(row(X)< col(X))
  X[ind2] <- x
  X[ind] <- X[ind]/sqrt(2)
  X[ind3] <- t(X)[ind3]
  return(X)
}


######################
######################
# TO-DO: use efficient eigen-decomposition due to dpr1 structure
eigenA.old <- function(w, N, r, tind, nP1, nq1){
  z <- rep(-1, N)
  z[tind] <- w-1
  z <- t(nP1) %*% z
  return(eigen(z %*% t(z) + diag(nq1)))
}

eigenA <- function(w, N, r, tind, nP1, nq1){
  z <- rep(-1, N)
  z[tind] <- w-1
  z <- t(nP1) %*% z
  return(dpr1_eigen(nq1, z, need.order=F)) # as nq1 is ordered
}

check.wbound <- function(w, lower, upper, tol=1e-8){
  as.logical((abs(lower-w) <= tol) + (abs(upper-w) <= tol))
}

#######################
#### main function ####
#######################

mineig <- function(w0, N, r, tind, nP1, nq1, lam2, lower=rep(1,length(w0)),
                   upper=rep(Inf, length(w0)), rho0=1e-2, tau0=1e-6, tol=1e-8,
                   max.iter=500, verbose=F){
  # step 0
  rho <- rho0
  tau <- tau0
  w <- w0
  m <- length(w)
  nn <- length(nq1)
  indices <- which(tind)
  z <- rep(-1,N)

  iter <- 1
  steps <- array(dim=max.iter)
  tts <- array(dim=max.iter)
  eA <- eigenA(w, N, r, tind, nP1, nq1)
  eA$values <- eA$values + lam2/N*sum(w^2)
  z[tind] <- w-1

  while (iter <= max.iter){
    if (verbose) cat("Iteration", iter, "\n")

    # step 1
    tt <- sum((eA$values[1] - eA$values) <= tau * max(1, abs(eA$values[1])))
    Q1 <- eA$vectors[, 1:tt, drop=F]
    tts[iter] <- tt

    ## linear programming
    obj <- c(1, rep(0,m))

    B <- nP1 %*% Q1

    tt2 <- (tt*(tt+1)/2)
    M <- matrix(0, nr=(tt2 + (nn-tt)), nc=(1+m))
    rhs <- array(dim=nrow(M))

    # (25), (26)
    M[1:tt2,1] <- svec(diag(tt))
    M[tt2+(1:(nn-tt)),1] <- 1

    for (k in (1:m)){
      kk <- indices[k]
      dA <- matrix(0, nr=N, nc=N)
      dA[kk,] <- z; dA[,kk] <- z; dA[kk,kk] <- 2*z[kk]
      M[1:tt2,1+k] <- -svec(t(B) %*% dA %*% B + 2*lam2/N * w[k] * diag(tt))
      PdAP <- t(nP1) %*% dA %*% nP1
      for (l in (1:(nn-tt))){
        M[tt2+l, 1+k] <- -as.vector(t(eA$vectors[, tt+l]) %*% PdAP %*%
                                    eA$vectors[, tt+l] + 2*lam2/N * w[k])
      }
    }
    rhs[1:tt2] <- svec(diag(eA$values[1:tt] - eA$values[1], nrow=tt, ncol=tt))
    rhs[tt2+(1:(nn-tt))] <- eA$values[(tt+1):nn] - eA$values[1]

    # (29), (27)
    # assume w is feasible, then (27) and (29) together: pmax(l-w, -rho) <= d <= pmin(u-w, rho)
    # use option: bounds
    bounds <- list(lower=list(ind=1:(m+1), val=c(-Inf, pmax(lower-w, -rho))),
                   upper=list(ind=1:(m+1), val=c(Inf, pmin(upper-w, rho))))

    # solve linear programming
    dirr <- c(rep("==", tt2), rep(">=", nn-tt))
    lpres <- Rglpk:::Rglpk_solve_LP(obj=obj, mat=M, dir=dirr, rhs=rhs, bounds=bounds,
                            control=list("verbose"=F))

    # step 2
    # check semi-positivity of U
    U <- svec.inv(lpres$auxiliary$dual[1:tt2])
    eU <- eigen(U)
    npd <- (eU$values[tt] < -1e-8) # check if the smallest eigenvalue negative
    if (npd) tau <- tau/2

    if (sqrt(sum(lpres$solution[-1]^2)) > tol){
      eA1 <- eigenA(w+lpres$solution[-1], N, r, tind, nP1, nq1)
      eA1$values <- eA1$values + lam2/N*sum((w+lpres$solution[-1])^2)
      if (eA1$values[1] >= eA$values[1]){
        # step 3
        if (verbose) cat("     step 3: rho reduced by 2\n")
        rho <- rho / 2   # update rho
        steps[iter] <- 3
      } else {
        # step 4
        if (verbose) cat("     step 4: update\n")
        psi <- (eA$values[1] - eA1$values[1])/(-lpres$solution[1])
        if (psi > 0.75) rho <- rho * 2
        if (psi < 0.25) rho <- rho / 2
        w <- w + lpres$solution[-1]
        z[tind] <- w-1
        eA <- eA1
        steps[iter] <- 4
      }
    } else {
      if (npd){
        # step 5 (eigenvalue splitting, Thm 7)
        if (verbose) cat("     step 5: eigenvalue splitting\n")
        eU$vectors[,tt]
        dind <- check.wbound(w, lower, upper)
        sol <- qr.solve(M[1:tt2, c(T, dind)], svec(eU$vectors[,tt] %*% t(eU$vectors[,tt])))
        d <- rep(0, m); d[dind] <- sol[-1]
        w <- w + d
        z[tind] <- w-1
        eA <- eigenA(w, N, r, tind, nP1, nq1)
        eA$values <- eA$values + lam2/N*sum(w^2)
        tau <- tau / 10
        steps[iter] <- 5
      } else {
        steps[iter] <- 0
        break # check if break correctly
      }
    }
    iter <- iter + 1
  }
  return(list(w=w, iter=iter, rho=rho, tau=tau, tt=tt, steps=steps, tts=tts))
}
