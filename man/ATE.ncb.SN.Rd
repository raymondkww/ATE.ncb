\name{ATE.ncb.SN}
\alias{ATE.ncb.SN}
\title{
  Kernel-based covariate balancing with lambda selection
}
\description{
  This function obtains the weights described in Wong and Chan (2018).
}
\usage{
ATE.ncb.SN(ind, K, lam1s, lam2s=1e-2*lam1s, lower=1, upper=Inf,
           thresh.ratio=1e-8, traceit=TRUE, w0=NULL, maxit=2000, maxit2=200,
           xtol_rel=1e-8, xtol_rel2=1e-4, method=2, check=FALSE, full=FALSE)
}
\arguments{
  \item{ind}{
    indicator vector of observation (T=1)
  }
  \item{K}{
    Gram matrix
  }
  \item{lam1s}{
    vector of lambda1
  }
  \item{lam2s}{
    vector of lambda2
  }
  \item{lower}{
    lower bound of weights
  }
  \item{upper}{
    upper bound of weights
  }
  \item{thresh.ratio}{
    threshold ratio for eigenvalue of K
  }
  \item{traceit}{
    print results or not
  }
  \item{w0}{
    initial value of weights
  }
  \item{maxit}{
    maximum number of iterations for BFGS
  }
  \item{maxit2}{
    maximum of iterations for SLP
  }
  \item{check}{
    check if max eigenvalue has multiplicity and, if so, apply SLP algorithm
  }
  \item{full}{
    return the full optimization results (reslist, reslist2) or not
  }
}
\references{R. K. W. Wong and K. C. G. Chan. (2018) "Kernel-based Covariate Functional Balancing for Observational Studies". Biometrika, 105(1), 199-213.}
\examples{
#######################
#### simulate data ####
#######################
set.seed(15)
n <- 200
Z <- matrix(rnorm(4*n),ncol=4,nrow=n)
prop <- 1 / (1 + exp(Z[,1] - 0.5 * Z[,2] + 0.25*Z[,3] + 0.1 * Z[,4]))
treat <- rbinom(n, 1, prop)
Y <- 200 + 10*treat+ (1.5*treat-0.5)*(27.4*Z[,1] + 13.7*Z[,2] +
                                      13.7*Z[,3] + 13.7*Z[,4]) + rnorm(n)
X <- cbind(exp(Z[,1])/2,Z[,2]/(1+exp(Z[,1])),
           (Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)
EY1X <- 200 + 10+ (1.5-0.5)*(27.4*Z[,1] + 13.7*Z[,2] +
                             13.7*Z[,3] + 13.7*Z[,4])
EY0X <- 200 + (-0.5)*(27.4*Z[,1] + 13.7*Z[,2] +
                      13.7*Z[,3] + 13.7*Z[,4])

w0 <- 1/prop*treat + 1/(1-prop)*(1-treat) # inverse propensity


mean(w0*treat*Y)-mean(w0*(1-treat)*Y) # ATE estimate based on inverse propensity (truth=10)


###########################################
##### kernel-based covariate balancing ####
###########################################

#### T=1 ####

# Sobolev kernel
Xstd <- transform.sob(X)$Xstd # standardize X to [0,1]^p
K <- getGram(Xstd) # get Gram matrix using Sobolev kernel

# design a grid for the tuning parameter
nlam <- 50
lams <- exp(seq(log(1e-8), log(1), len=nlam))

# compute weights for T=1
fit1 <- ATE.ncb.SN(treat, K, lam1s=lams)
if (sum(fit1$warns)) cat("lambda bound warning!\n")


#### T=0 ####

# compute weights for T=0
fit0 <- ATE.ncb.SN(1-treat, K, lam1s=lams)
if (sum(fit0$warns)) cat("lambda bound warning!\n")


#### ATE ####
mean(fit1$w*Y - fit0$w*Y) # ATE estimate based on kernel-based estimation (truth=10)

}
