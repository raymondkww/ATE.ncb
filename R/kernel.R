#################
#### kernels ####
#################

#########################################
#########################################

K.gauss <- function(x,y,sigma=1){
  return(exp(-sum((x-y)^2) / sigma^2))
}

#########################################
#########################################

# standardize L2-norm=1
K.gauss2 <- function(x,y,sigma=1){
  return(dnorm( sqrt(sum((x-y)^2)), mean=0, sd=sigma))
}

#########################################
#########################################

K.laplace <- function(x,y,alpha){
  return(exp(-alpha*sum(abs(x-y))))
}


################################
##### Sobolev RK function ######
################################
# not optimized

#########################################
#########################################

k1 <- function(t){
  return(t-.5)
}

#########################################
#########################################

k2 <- function(t){
  return( (k1(t)^2-1/12)/2 )
}

#########################################
#########################################

k4 <- function(t){
  return( (k1(t)^4-k1(t)^2/2+7/240)/24 )
}

#########################################
#########################################

K.sob <- function(s,t){
  ans <- 1 + k1(s)*k1(t) + k2(s)*k2(t) - k4(abs(s-t))
  return(ans)
}

#########################################
#########################################

K.sob.prod <- function(x, y){
  p <- length(x)
  out <- 1
  for (j in (1:p)){
    out <- out*K.sob(x[j], y[j])
  }
  return(out)
}

#########################################
#########################################

standard <- function(x){
  return((x-mean(x))/sd(x))
}

#########################################
#########################################

transform.sob <- function(X){
  Xlim <- apply(X, 2, range)
  Xstd <- matrix(nr=nrow(X), nc=ncol(X))
  for (i in (1:ncol(X))){
    Xstd[,i] <- (X[,i]-Xlim[1,i])/diff(Xlim[,i])
  }
  return(list(Xstd=Xstd, Xlim=Xlim))
}


##################################################
#### general function to produce Gram matrix #####
##################################################

getGram <- function(X, standardize=T, ker="sob", sigma=1, alpha=1){
  # note that sobolev kernel is not compatible with standardize=T, and it requires X to be
  # transformed to [0,1]^p
  if (ker=="sob") standardize <- F
  if (standardize) X <- apply(X, 2, standard)
  N <- nrow(X)
  K <- matrix(nr=N, nc=N)
  for (i in (1:N)){
    for (j in (1:N)){
      if (ker=="gauss"){
        K[i,j] <- K.gauss(X[i,], X[j,], sigma=sigma)
      } else if (ker=="gauss2"){
        K[i,j] <- K.gauss2(X[i,], X[j,], sigma=sigma)
      } else if (ker=="laplace"){
        K[i,j] <- K.laplace(X[i,], X[j,], alpha=alpha)
      } else if (ker=="sob"){
        K[i,j] <- K.sob.prod(X[i,], X[j,])
      }
    }
  }
  return(K)
}


