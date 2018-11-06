#################
#### ATE.ncb ####
#################

#########################################
#########################################

eval_obj_grad <- function(w, N, r, tind, nP1, nq1, lam2){
  z <- rep(-1, N)
  z[tind] <- w-1
  V <- sum(w^2)/N
  z <- t(nP1) %*% z
  dpr1 <- dpr1_eigval(r, nq1, z, need.order=F, vec=T)
  v <- dpr1[[4]]
  #v <- eigen(z %*% t(z) + diag(nq1))$vectors[,1]
  return(list("objective"=dpr1[[2]] + lam2*V,
              "gradient"=as.vector(2 * as.vector(t(v) %*% z) * nP1[tind,] %*%
                                   v) + 2*lam2*(w)/N)) # modified
}

#########################################
#########################################

#' ATE.ncb core function
#'
#' @param ind indicator vector of observation (T=1)
#' @param K Gram matrix
#' @param lam1s vector of lambda1
#' @param lam2s vector of lambda2
#' @param lower lower bound of weights
#' @param upper upper bound of weights
#' @param thresh.ratio threshold ratio for eigenvalue of K
#' @param traceit trace it or not
#' @param w0 initial value of weights
#' @param maxit maximum number of iterations for BFGS
#' @param maxit2 maximum of iterations for SLP
#' @param check check if max eigenvalue has multiplicity and, if so, apply SLP algorithm
#' @param full return the full optimization results (reslist, reslist2)?


ATE.ncb.core <- function(ind, K, lam1s, lam2s=1e-2*lam1s, lower=1, upper=Inf,
                         thresh.ratio=1e-8, traceit=TRUE, w0=NULL, maxit=2000,
                         maxit2=200, xtol_rel=1e-8, xtol_rel2=1e-4,
                         check=FALSE, full=FALSE){ 
  N <- length(ind)
  
  # construct a vector
  tind <- as.logical(ind)
  n1 <- sum(tind)
  n2 <- sum(!tind)

  # lower and upper bound
  if (length(lower)==1) lower <- rep(lower, n1)
  if (length(upper)==1) upper <- rep(upper, n1)

  # construct P1 and q1
  e <- eigen(K)
  thresh <- e$values[1]*thresh.ratio
  eind <- (e$values>=thresh)
  r <- sum(eind)
  
  if (traceit)
    cat("Number of components used: r = ", r, "with threshod.ratio =", thresh.ratio, "\n")

  if (is.null(w0)){
    w0 <- lower+1
  }

  nlams <- length(lam1s)
  lam1s <- sort(lam1s) # lams is sorted
  outlist <- list()
  ws <- matrix(0, nr=nlams, nc=N)
  SNs <- array(dim=nlams)
  unorms <- array(dim=nlams)
  fittedus <- array(dim=c(nlams, nc=N))
  rmseTo1s <- array(dim=nlams)
  outlist <- NULL
  outlist2 <- NULL
  alg <- rep(1,nlams)
  obj1 <- array(dim=nlams)
  obj2 <- array(dim=nlams)

  if (length(lam2s)==1){
    lam2s <- rep(lam2s, nlams)
  }
  
  for (i in (1:nlams)){
    nP1 <- e$vectors[,eind]/sqrt(N) # P1/sqrt(N)
    nq1 <- -lam1s[i]*N/e$values[eind] # -lam1s[i]*N1*q1
    # reorder
    oo <- order(nq1)
    nq1 <- nq1[oo]
    nP1 <- nP1[,oo]

    if (traceit) cat("####", i, ":", "\tlam1 =", lam1s[i], "\tlam2 =", lam2s[i], "\n")

    res <- nloptr:::nloptr(x0=w0, eval_f=eval_obj_grad, lb=lower, ub=upper,
                           N=N, r=r, tind=tind, nP1=nP1, nq1=nq1, lam2=lam2s[i],
                           opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel=xtol_rel,
                                     print_level=0, maxeval=maxit, check_derivatives=F))

    if (full)  outlist[[i]] <- res
    obj1[i] <- res$objective
    ws[i,tind] <- res$sol
    temp <- t(nP1) %*% (ws[i,]-1)
    ee <- eigen(temp %*% t(temp) + diag(nq1))

    ## checking if the largest eigenvalue has multiplicity 1
    if (check){
      if (abs(ee$values[1]-ee$values[2])/abs(ee$values[1]) < 1e-6){
        cat("multiplicity issue of largest eigenvalues\n")
        res2 <- mineig(w0=res$sol, N=N, r=r, tind=tind, nP1=nP1, nq1=nq1, lam2=lam2s[i],
                      lower=lower, upper=upper, rho0=1e-3, tau0=1e-6,
                      tol=xtol_rel2, max.iter=maxit2, verbose=T)
        if (full) outlist2[[i]] <- res2
        obj2[i] <- eval_obj_grad(res2$w, N=N, r=r, tind=tind, nP1=nP1, nq1=nq1, lam2=lam2s[i])$objective

        # udpate w
        ws[i,tind] <- res2$w
        temp <- t(nP1) %*% (ws[i,]-1)
        ee <- eigen(temp %*% t(temp) + diag(nq1))
        alg[i] <- 2
      }
    }

    # compute SN
    SNs[i] <- (sum(temp*ee$vectors[,1]))^2
    unorms[i] <- sqrt(sum(nq1/(-lam1s[i]) * ee$vectors[,1]^2))
    fittedus[i,] <- as.vector(nP1 %*% ee$vectors[,1] * (N))
    rmseTo1s[i] <- sqrt(mean((ws[i,]-1)^2))

    w0 <- res$sol
  }

  return(list(outlist=outlist, outlist2=outlist2, ws=ws, SNs=SNs, unorms=unorms,
              fittedus=fittedus, rmseTo1s=rmseTo1s, lam1s=lam1s, lam2s=lam2s, alg=alg, obj1=obj1, obj2=obj2,
              check=check))
}


#########################################
#########################################

#' ATE.ncb with lambda selection 
#'
#' @param ind indicator vector of observation (T=1)
#' @param K Gram matrix
#' @param lam1s vector of lambda1
#' @param lam2s vector of lambda2
#' @param lower lower bound of weights
#' @param upper upper bound of weights
#' @param thresh.ratio threshold ratio for eigenvalue of K
#' @param traceit trace it or not
#' @param w0 initial value of weights
#' @param maxit maximum number of iterations for BFGS
#' @param maxit2 maximum of iterations for SLP
#' @param check check if max eigenvalue has multiplicity and, if so, apply SLP algorithm
#' @param full return the full optimization results (reslist, reslist2)?


ATE.ncb.SN <- function(ind, K, lam1s, lam2s=1e-2*lam1s, lower=1, upper=Inf,
                       thresh.ratio=1e-8, traceit=TRUE, w0=NULL, maxit=2000,
                       maxit2=200, xtol_rel=1e-8, xtol_rel2=1e-4, method=2,
                       check=FALSE, full=FALSE){
  ores <- ATE.ncb.core(ind=ind, K=K, lam1s, lam2s=lam2s, lower=lower, upper=upper,
                       thresh.ratio=thresh.ratio, traceit=traceit, w0=w0,
                       maxit=maxit, maxit2=maxit2, xtol_rel=xtol_rel, xtol_rel2=xtol_rel2,
                       check=check, full=full)

  # find which SN is the smallest
  if (method==1){
    ind <- which.min(ores$SNs)
  } else if (method==2){
    ind <- which(diff(ores$SNs)/diff(ores$lam1s) > (-1e-6))[1]
  }
  return(list(w=ores$ws[ind,], ind=ind, warns=c(ind==1, ind==length(lam1s)), ores=ores))
}
