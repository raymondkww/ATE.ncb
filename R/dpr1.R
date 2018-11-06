# compute the i-th eigenvalue, with i=length(d) being the maximum eigenvalue,
# of the matrix diag(d) + outer(z, z)
# the algorithm requires both d and z are ordered according to d
# vec=T: egeinvector will be returned as res[[4]]

dpr1_eigval <- function(i, d, z, need.order=T, vec=F){
  if (need.order){
    o <- order(d)
    d <- d[o]
    z <- z[o]
  }
  rho <- sum(z^2)
  z <- z/sqrt(rho)
  res <- .Call("dpr1_eigval", as.integer(i), d, z, as.double(rho), PACKAGE="ATE.ncb")
  if (vec){
    vv <- 1/res[[1]]*z
    res[[4]] <- vv/sqrt(sum(vv^2))
  }
  return(res)
}


# compute the full eigendecomposition for matrix diag(d) + outer(z, z)
dpr1_eigen <- function(d, z, Kstart=1, Kstop=length(d), need.order=T, order.eigen=T){
  if (need.order){
    o <- order(d)
    d <- d[o]
    z <- z[o]
  }
  rho <- sum(z^2)
  z <- z/sqrt(rho)
  res <- .Call("dpr1_eigen", as.integer(Kstart), as.integer(Kstop), d, z, as.double(rho), PACKAGE="ATE.ncb")
  if (order.eigen){
    nn <- ncol(res[[1]])
    res[[1]] <- res[[1]][,nn:1]
    res[[2]] <- res[[2]][nn:1]
  }
  names(res) <- c("vectors", "values", "info")
  return(res)
}

