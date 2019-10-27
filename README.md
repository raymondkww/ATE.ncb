# ATE.ncb: Kernel-based covariate functional balancing for observational studies

This package implements the covariate balancing method proposed in Wong and Chan (2018).

## Installation
This package can be installed via function `install_github` in R package `devtools`:

``` r
install.packages("devtools")
devtools::install_github("raymondkww/ATE.ncb")

```

## Example
```r
library("ATE.ncb")

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
```


## References
* R. K. W. Wong and K. C. G. Chan. (2018) "Kernel-based Covariate Functional Balancing for Observational Studies". Biometrika, 105(1), 199-213. [\[link\]](https://doi.org/10.1093/biomet/asx069)
