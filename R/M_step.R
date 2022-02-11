## M_step.R

## Update pi.list, LLc, LLint, long.df, this.glm

M_step <- function(family, formula, response,
                   nf4c, nlev, nclus,
                   f4c.colno, cf.colno,
                   indcf.list, pp.list,
                   long.df, verbose){

    if (verbose > 1) cat("Starting M step\n")

    ## Update pi.list using updated pp.list:
    pi.list <- vector("list", nf4c)
    for (cf in 1:nf4c)
        pi.list[[cf]] <- apply(pp.list[[cf]],2,mean)

    names(pi.list[[1]]) <- paste("pi",1:nclus[1],sep="")
    if (nf4c == 2)
        names(pi.list[[2]]) <- paste("kappa",1:nclus[2],sep="")
    
    if (verbose > 1){
        cat("Updated pi.list\n")
        cat("[[1]]\n")
        print(pi.list[[1]])
        if (nf4c > 1){
            cat("[[2]]\n")
            print(pi.list[[2]])
        }
    }
        
    ## Put (or update) long pi, kappa columns in long.df:
    pi.longv <- c(indcf.list[[1]] %*% pi.list[[1]])
    long.df$pi <- pi.longv
    if (nf4c == 2){
        kappa.longv <- c(indcf.list[[2]] %*% pi.list[[2]])
        long.df$kappa <- kappa.longv
    }
        
    if (verbose > 1) print(pi.longv)
    
    ## Construct pp1 (pp2) and weights columns using values from pp.list:

    Nlong <- nrow(long.df)
    long.df$wts <- rep(1,Nlong)
    if ((family == "binomial")&(is.vector(response))&!is.null(long.df$weights))
        long.df$wts <- long.df$weights
    ## Already have prior weights from std.glm, either 1 or ntrials
    wts <- long.df$wts
    pp1 <- rep(1,Nlong)
    pp2 <- rep(1,Nlong)
    for (nn in 1:Nlong){
        
        ic <- as.numeric(as.character(long.df[nn,c(f4c.colno[1],
                                                   cf.colno[1])]))
        pp1[nn] <- pp.list[[1]][ic[1],ic[2]]
        wts[nn] <- wts[nn] * pp1[nn]
        if (nf4c > 1){
            ic <- as.numeric(as.character(long.df[nn,c(f4c.colno[2],
                                                       cf.colno[2])]))
            pp2[nn] <- pp.list[[2]][ic[1],ic[2]]
            wts[nn] <- wts[nn] * pp2[nn]
        }
    }

    long.df$pp1 <- pp1
    if (nf4c > 1) long.df$pp2 <- pp2

    ## Weights (and prior weights from binomial)
    ## Use matrix input of response 
    delt <- 0.000000001
    wts[wts<=delt] <- delt ## glm won't allow wts of zero STILL TO FIX   
    long.df$wts <- wts

    ## If binomial data with matrix input, glm will construct prior weights.
    ## Just make wts column the pp1 values, glm will multiply by ntrials.

    if (is.matrix(response)){
        long.df$wts <- long.df$pp1
        if (verbose > 1) print(long.df$wts)
    }
    
    ## Fit the glm to the long data:
    
    this.glm <- glm(formula=formula,weights=wts,
                    data=long.df,family=family)

    ## Ancilliary parameters:        
    
    if (family=="gaussian"){
	## Constant variance, scalar:
        var.s <- this.glm$deviance/this.glm$df.residual
	long.df$sigsq <- rep(var.s,Nlong)
    }

    ## if (family == "nbin") Variance inflation, COMPLETE LATER        

    ## Find log likelihood with complete data, LLc, using glm().
    ## Deviance more reliable than logLik over different distributions
    ## logLik term 1:
    LLc <- -this.glm$deviance/2
    ## Adjust for multinomials:
    for (cf in 1:nf4c) LLc <- LLc +
        (nclus[cf]>1)*(nclus[cf]<nlev[cf])* ## Ensure non-trivial
            nlev[cf]*sum(xlogx.fn(pi.list[[cf]]))
    ## This is LLc = log likelihood under complete knowledge

    ## Update mu = conditional fitted value in raw data scale,
    long.df$mu <- fitted(this.glm)  ## OK, assuming no discarded rows

    if (verbose > 1){
        print("long.df$mu")
        print(long.df$mu)
        print("names of long.df")
        print(names(long.df))
        print(long.df[1:3,])
    }
    
    ## Update pdf(y,theta) and logpdf(y,theta) columns:        
    
    if (family=="binomial"){        
        long.df$logpdf <- binomial.logpdf(long.df$mu, long.df$ntrials, long.df$nsucc)
        long.df$pdf <- exp(long.df$logpdf)        
    }        
    if (family=="poisson"){        
        long.df$logpdf <- poisson.logpdf(long.df$mu,long.df$y)        
        long.df$pdf <- exp(long.df$logpdf)        
    }        
    if (family=="gaussian"){        
        long.df$logpdf <- gaussian.logpdf(long.df$mu,long.df$sigsq,long.df$y)        
        long.df$pdf <- exp(long.df$logpdf)        
    }        
    
    if (verbose > 1){
        print("long.df$logpdf")
        print(long.df$logpdf)
        print("long.df$pdf")
        print(long.df$pdf)
    }
    
    ## Update pi*pdf and pp1*pdf columns:
    long.df$pi.pdf <- long.df$pi * long.df$pdf
    long.df$pp1.pdf <- long.df$pp1 * long.df$pdf        

    if (verbose > 1){
        print("long.df$pi.pdf")
        print(long.df$pi.pdf)
        print("long.df$pp1.pdf")
        print(long.df$pp1.pdf)
    }

    ## If biclustering, update kappa*pdf and pp2*pdf columns:
    if (nf4c == 2){        
	long.df$kappa.pdf <- long.df$kappa * long.df$pdf        
	long.df$pp2.pdf <- long.df$pp2 * long.df$pdf        
    }        
    
    ## Find integrated (marginal, incomplete) LL:
    
    LLint <- LLc

    ## Add entropy term(s):

    LLint <- LLint - sum(xlogx.fn(c(pp.list[[1]])))

    if (length(nf4c) == 2)
        LLint <- LLint - sum(xlogx.fn(c(pp.list[[2]])))

    return(list(pi.list = pi.list, pp.list = pp.list,
                LLc = LLc, LLint = LLint,
                long.df = long.df, this.glm = this.glm))
}
