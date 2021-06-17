## ---------
## somefun.R
## ---------

## Some functions for clustglm.


## Function names:
## ---------------

## replacedefaults

## xlogx.fn
## log.fn
## logit.fn
## expit.fn
## id.fn

## maxmiss.fn

## binomial.logpdf
## poisson.logpdf
## NB.logpdf
## gaussian.logpdf

## AIC.fn
## AICc.fn
## BIC.fn

## run.clustglm
## mat2df
## longdf2mat

## ICtable.fn


#####################################################

## replacedefaults
## ---------------

replacedefaults <- function (default, user) replace(default, names(user), user)



#####################################################

## xlogx function:
## ---------------

## Function for x.log(x) with value zero if x=0 
## (or if computer puts x just below zero). 
## Argument may be a scalar, vector or matrix.

xlogx.fn <- function(x){
    xlx <- x
    xlx[xlx <= 0] <- 0
    xlx[xlx > 0]  <- x[x>0] * log(x[x>0])
    return(xlx)
}


## log function:
## -------------

## If p is really close to 0, set log to some min.
## This min is around -20.

log.fn <- function(x){
    ## Replace values if too near 0:
    x[x<0.000000001] <- 0.000000001
    ## Output object:
    return(log(x))
}


## logit function:
## ---------------

## If p is really close to 0 or 1, set logit to some max or min.

logit.fn <- function(x){
    ## Replace values if too near 0 or 1:
    x[x<0.000000001] <- 0.000000001
    x[x>0.999999999] <- 0.999999999
    ## Output vector:
    return(log(x/(1-x)))
}

## expit function:
## ---------------

## Inverse of logit function, acts on a logit vector.

expit.fn <- function(x)
    return(1/(1+exp(-x)))


## identity function:
## ------------------

id.fn <- function(x)
    return(x)


## maxmiss function. deals with missing data points in E step, stops warnings:
## ---------------------------------------------------------------------------

maxmiss.fn <- function(x){
    nmiss <- length(x[is.na(x)])
    if (nmiss<length(x)) y <- max(x,na.rm=T)
    if (nmiss==length(x)) y <- -Inf
    return(y)
}
    

## randpart.fn, generates a random partition of n objects in c classes
## with non-empty classes

randpart.fn <- function(n,c){
    pp <- matrix(0,n,c)
    pp[1:c,] <- diag(rep(1,c))  # Ensures non-empty classes
    pp[(c+1):n,1] <- 1
    for (i in (c+1):n)
        pp[i,] <- sample(pp[i,])
    pp <- pp[sample(1:n),]
    pp
}
                                        #

## logpdf functions:
## -----------------

binomial.logpdf <- function(theta.v, nt.v, y.v){
    lgamma(nt.v+1) - lgamma(y.v+1) - lgamma(nt.v-y.v+1) +
        y.v*log.fn(theta.v) + (nt.v-y.v)*log.fn(1-theta.v)
}

poisson.logpdf <- function(mu.v,y.v){
    -mu.v + y.v * log.fn(mu.v) - lgamma(y.v+1)
}

NB.logpdf <- function(mu.v,k.v,y.v){
    lpdf <- lgamma(y.v + k.v) - lgamma(y.v+1) -
        lgamma(k.v) + y.v*log.fn(mu.v) +
            k.v*log.fn(k.v) -
                (y.v+k.v)*log.fn(mu.v+k.v)
    return(lpdf)
}

gaussian.logpdf <- function(mu.v,var.v,y.v){
    -log(2*pi*var.v) -0.5*(y.v-mu.v)^2/var.v
}


