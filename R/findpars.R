## findpars.R

## Find parameters using tapply and weighted sum-to-zero
## constraints on the linear predictors from the fitted model.
## The fitted model may have at most one interaction term.
## Each type of parameter allows for previous ones in model,
## i.e. a sequential analysis of deviance (Type I rather than
## Type III SS in SAS).

## The glm() used inside M.step of clustglm() uses corner-point
## parameterisation, zero for first level of each factor.
## This findpars function uses weighted sum-to-zero constraints,
## providing interpretable parameters for profile plots,
## ordination plots and biplots. This choice of constraints and
## the resultant plots are most suited to balanced data, although
##  for covariates and unbalanced data.

## Standard errors are not yet implemented.


findpars <- function(model){
    this.model <- model
    long.df <- this.model$data
    N <- nrow(long.df)    ## 

    my.form <- this.model$formula  ## e.g. Y ~ A+B+C:D, or Y ~ C:D
    fterms <- terms(my.form)       ## Tried keep.order = TRUE
    ftermlabels <- attr(fterms, "term.labels")
    ftermorder <- attr(fterms, "order")
    ftermfactors <- rownames(attr(fterms, "factors"))
    fterms.split <- strsplit(ftermlabels, ":")
    ## R may have reordered the interaction components;
    ## using keep.order didn't stop it

    ## Check at most one two-way interaction:
    if (sum(ftermorder == 2) > 1)
        stop("This function requires a model with at most one two-way interaction")

    ## Set up list of main effects parameters.
    nME <- sum(ftermorder == 1)   ## Number of main effects; may be zero
    if (nME > 0){
        ME.list <- vector("list", nME)
        names(ME.list) <- paste("ME", 1:nME, sep = "")
        for (me in 1:nME){
            this.name <- fterms.split[[me]]
            names(ME.list)[me] <- this.name
            this.factor <- long.df[,(1:ncol(long.df))[this.name == names(long.df)]]
            this.nlevels <- length(levels(this.factor))
            ME.list[[me]] <- list(factor = this.factor,
                                  nlevels = this.nlevels,
                                  pars = rep(NA,this.nlevels),
                                  pars.long = rep(NA,N))
        }    

    }

    ## If there is an interaction term:
    if (sum(ftermorder == 2) == 1){
        INT1 <- fterms.split[[nME+1]][1]
        INT2 <- fterms.split[[nME+1]][2]
        ## Were the interaction labels reversed?
        if (paste(INT2,":",INT1, sep = "") %in%
            strsplit(as.character(my.form)," ")[[3]]){
            temp <- INT1
            INT1 <- INT2
            INT2 <- temp
        }
        ## Set up interaction factors and levels:
        INT1.f <- long.df[,(1:ncol(long.df))[INT1 == names(long.df)]]
        INT1.nlevels <- length(levels(INT1.f))
        INT2.f <- long.df[,(1:ncol(long.df))[INT2 == names(long.df)]]
        INT2.nlevels <- length(levels(INT2.f))
    }

    ## Set up output list
    pars.out <- vector("list", 1 + length(ftermorder))
    if (length(ftermorder) == 4) if (all(ftermorder == c(1,1,1,2)))
        names(pars.out) <- c("nu","alpha","beta","psi","gamma")
    if (length(ftermorder) == 3) if (all(ftermorder == c(1,1,2)))
        names(pars.out) <- c("nu","alpha","beta","gamma")
    if (length(ftermorder) == 2) if (all(ftermorder == c(1,2)))
        names(pars.out) <- c("nu","alpha","gamma")
    if (length(ftermorder) == 1) if (all(ftermorder == c(2)))
        names(pars.out) <- c("nu","gamma")
    if (length(ftermorder) == 3) if (all(ftermorder == c(1,1,1)))
        names(pars.out) <- c("nu","alpha","beta","psi")
    if (length(ftermorder) == 2) if (all(ftermorder == c(1,1)))
        names(pars.out) <- c("nu","alpha","beta")
    if (length(ftermorder) == 1) if (all(ftermorder == c(1)))
        names(pars.out) <- c("nu","alpha")
    if (length(ftermorder) == 0)
        names(pars.out) <- c("nu")
    if (length(ftermorder) > 0){
        if (sum(ftermorder == 2) == 1)
            names(pars.out) <- c("Overall mean", names(ME.list),
                                 paste(INT1,":",INT2, sep = ""))
        if (sum(ftermorder == 2) == 0)
            names(pars.out) <- c("Overall mean", names(ME.list))
    }

    ## Use tapply() to find parameters:
    ## If unclustered model, assign weights 1:
    if (is.null(long.df$wts)){
        Nobs <- nrow(long.df)
	long.df$wts <- rep(1,Nobs)
    }
    lp <- this.model$linear.predictors
    wlp <- long.df$wts*lp
    wts <- long.df$wts

    ## nu, overall weighted mean
    nu <- sum(wlp)/sum(wts)
    pars.out[[1]] <- nu
    
    ## Main effects parameters
    if (nME > 0){
        for (me in 1:nME){
            par.v <- tapply(wlp, ME.list[[me]]$factor, sum)/
                tapply(wts, ME.list[[me]]$factor, sum) - nu
            names(par.v) <- levels(ME.list[[me]]$factor)
            long.v <- rep(NA,N)
            for (i in 1:length(par.v))
                long.v[names(par.v[i]) ==
                      ME.list[[me]]$factor] <- par.v[i]
            ME.list[[me]]$pars <- par.v
            ME.list[[me]]$pars.long <- long.v
            ## Save for output:
            pars.out[[me+1]] <- par.v
        }
    }

    ## Interaction parameters
    ## Calculate the gamma matrix for this model:
    if (sum(ftermorder == 2) == 1){
        ## Starting term:
        num <- tapply(wlp, list(INT1.f,INT2.f), sum)
        ## Subtract main effects terms:
        if (nME > 0) for (me in 1:nME){
            num <- num - tapply(wts*ME.list[[me]]$pars.long,
                                list(INT1.f,INT2.f), sum)
        }
        den <- tapply(wts, list(INT1.f,INT2.f), sum)
        gam.m <- num/den - nu
        pars.out[[1+nME+1]] <- gam.m
    }


    return(pars.out)
}


 
