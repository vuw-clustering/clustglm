## clustglm.R
## -----------
## Version 1.4
## clustglm function only
## See somefun.R for other functions, including replacedefaults

######################################################



## Model fitting function:
## -----------------------

clustglm <- function(formula, 
                     family = c("gaussian", "binomial", "poisson"), 
                     data,
                     weights = NULL,  ## Prior weights, e.g. ntrials for binomial
                     offset,
                     control = list(...),
                     
                     ## Clustered model specification:
                     fact4clust = NULL,       # Factors to be clustered
                     nclust = NULL,            # Numbers of clusters
                     clustfactnames = NULL,   # New names for clustered factors
                     
                     ## Options for finding starting points:
                     start.control = list(randstarts = 0, selfstarts = 0, alloc = NULL),

                     ## Control of EM algorithm:
                     EM.control = list(maxEMcycles = 100, EMstoppingpar = 1e-04),

                     ## Save long data frame:
                     save.long = TRUE,
                     
                     ## Save sequence of estimates:
                     save.ests = FALSE,

                     ## Save quantile residuals (continuous family, e.g. Gaussian):
                     save.Qres = FALSE,

                     ## Save randomised quantile residuals (discrete family, e.g. Poisson):
                     save.RQres = FALSE,

                     ## Save expected quantile residuals (discrete family, e.g. Poisson):
                     save.EQres = FALSE,

                     verbose = 0,    ## Can change to 1 or 2 for medium or high verbosity
                     ...
                     ) {

    ## Replace defaults with user-provided values:
    ## -------------------------------------------
    
    default.start.control = as.list(args(clustglm))$start.control
    default.EM.control = as.list(args(clustglm))$EM.control
    start.control <- replacedefaults(default.start.control, start.control)
    EM.control <- replacedefaults(default.EM.control, EM.control)
    
    ## Sort out weights for binomial:
    ## ------------------------------

    if ((family == "binomial")&(!(is.null(weights)))){
        wtcolno <- (1:ncol(data))[weights == names(data)]
	## or    wtcolno <- which(names(data) == weights)
        data$weights <- data[,wtcolno]
    }

    ## Check names of factors and response variable(s)
    ## -----------------------------------------------
    
    fo <- formula

    response <- NULL  ## Default for Binomial with matrix data entry
    ## If single response variable, find column number, check OK:
    if (class(attr(terms(fo),"variables")[[2]]) == "name"){
        response.name <- attr(terms(fo),"variables")[[2]]
        if (as.character(response.name) %in% names(data)){
            response.colno <- (1:ncol(data))[names(data) == as.character(response.name)]
            if (length(response.colno) > 1)
                stop("More than one column in the data set has the response label")
            ## Check the single response is numeric:
            response <- data[,response.colno]
            if ((length(response.colno) == 1)&(!(is.numeric(response))))
  	        stop("Response variable must be numeric")
	    ## If OK, keep this response.
        }
    }

    ## Do more checks:
    if (family == "gaussian"){
        if ((save.EQres == TRUE)|(save.RQres == TRUE)){
            cat("EQres and RQres are for discrete data.\n")
            stop("Use save.Qres for continuous data (e.g. Gaussian).\n")
            }
    }
            
    if (family == "poisson"){
        if (!((all(as.integer(response) == response)) & (min(response) >= 0)))
            stop("All Poisson responses must be non-negative integers")
	if (save.Qres == TRUE){
	    cat("Qres is for continuous distributions, Poisson is discrete.\n")
	    stop("Use save.RQres and/or save.EQres instead.")
	}
    }
    if (family == "binomial"){
	if (save.Qres == TRUE){
	    cat("Qres is for continuous distributions, Binomial is discrete.\n")
	    stop("Use save.RQres and/or save.EQres instead.")
	}
        ## Different types of data entry:
        ## (i) Response data binary 0/1:
        if (is.numeric(response)){
            if (all(response %in% c(0,1))){
                data$nsucc <- response
                data$ntrials <- rep(1, nrow(data))
                data$nfail <- 1 - data$nsucc
                data$psucc <- data$nsucc / data$ntrials
                if (verbose > 0)
                    cat("NOTE: All responses are 0 or 1; we assume binary data from single trials\n")
            }
        }
        ## (ii) Response data on interval [0,1], prob. success.
        ## Check if outside interval.
        if (is.numeric(response)){
            if ((max(response)>1) | (min(response) < 0)){
                stop("All responses should be probabilities, in the interval [0,1]")
            }
	}
        ## (iii) Response data is a call: cbind(nsucc, nfail).
        if (class(attr(terms(fo),"variables")[[2]]) == "call"){
            response.call <- attr(terms(fo),"variables")[[2]]
	    attach(data)
	    response <- eval(response.call)  ## Matrix or ratio
	    detach()
            ## Two types of calls, matrix or ratio.
            if (is.matrix(response)){
                data$nsucc <- round(response[,1])
                data$nfail <- round(response[,2])
                data$ntrials <- apply(response,1,sum)
                data$weights <- data$ntrials
                data$psucc <- data$nsucc / data$ntrials
            }
	}
    }
    if (verbose > 1){
        cat("Names of the data frame \n")
        cat(names(data), "\n")
    }

    ## -------------------------------------------------------------
    ## Run a preliminary glm, using unclustered terms in the formula
    ## -------------------------------------------------------------

    ## Create a simplified formula, removing clustered terms, then fit fixed.glm:
    ## Read in and check the formula:
    glmformula <- formula
    ## drop any terms from glmformula that contain variables not in data
    fterms <- terms(formula)
    ftermlabels <- attr(fterms, "term.labels")  
    if (length(ftermlabels)>0) {
        OKterm.fn <- function(x) !all(x %in% names(data))
        dropi <- which(sapply(strsplit(ftermlabels, ":"), OKterm.fn))
        if ((length(dropi) > 0)&(length(dropi)<length(ftermlabels))) {
            glmterms <- drop.terms(fterms, dropi, keep.response = TRUE)
            glmformula <- reformulate(attr(glmterms,"term.labels"), glmterms[[2]])
        }
        if ((length(dropi) > 0)&(length(dropi) == length(ftermlabels))) {
	    ## Losing all terms, need intercept for null model
	    glmformula <- response ~ 1
        }
    }

    ## Fit the preliminary glm (fixed effects only):
    if (is.null(weights)) fixed.glm <- glm(formula = glmformula, family = family,
                                         data = data)
    if (!(is.null(weights))) fixed.glm <- glm(formula = glmformula, family = family,
                                         data = data, weights = weights)
    
    ## Save the response as y, needed in M step.
    data$y <- fixed.glm$y
    
    if (verbose > 0) cat("Preliminary glm fitted\n")

    ## Set up default prior weights:
    Nobs <- nrow(data)
    data$wts <- rep(1,Nobs)

    ## Store binomial matrix input under new names:
    if ((family == "binomial") & (is.matrix(response))){
        data$psucc <- fixed.glm$y
        data$ntrials <- fixed.glm$weights
        data$nsucc <- round(data$psucc * data$ntrials)
        data$nfail <- data$ntrials - data$nsucc
    }
    ## Have checked, this works

    ## Store fixed.glm (30 items) then modify for correct LL:
    model.out <- fixed.glm
    model.out$npar <- fixed.glm$rank
    model.out$LL <- -fixed.glm$deviance/2
    model.out$AIC <- fixed.glm$deviance + 2*fixed.glm$rank
    model.out$BIC <- fixed.glm$deviance + fixed.glm$rank*log(sum(data$wts))

    if ((save.Qres == TRUE)&(family == "gaussian")){
        ## Continuous distribution, use PIT (prob. integral transform)
        Qres <- pnorm(q = fixed.glm$y, mean = fixed.glm$fit, sd = fixed.glm$res)
        ## STILL TO TEST THIS FOR GAUSSIAN
        model.out$Qres <- Qres
    }
    
    ## Calculate RQres, EQres if requested, or if selfstarts requested:
    if (((save.RQres == TRUE)|(save.EQres == TRUE))|(start.control$selfstarts > 0)){
        ## Should be discrete distribution, negbin still to do
        ## Find lower and upper limits of risers
        lower <- rep(NA,Nobs)
        upper <- rep(NA,Nobs)
        mu <- model.out$fitted.values
        ## Assuming Poisson family
        if (family == "poisson"){
            y <- model.out$y
            MM <- max(y)
            for (nn in 1:Nobs){
                yy <- y[nn]
                mun <- mu[nn]
                dpois1 <- dpois(0:MM, lambda = mun)
                cumdpois1 <- c(0,cumsum(dpois1),1)
                lower[nn] <- cumdpois1[yy+1]
                upper[nn] <- cumdpois1[yy+2]
            }
        }
        ## Assuming binomial family:
        if (family == "binomial"){
            nsucc <- model.out$data$nsucc
            ntrials <- model.out$data$ntrials
            for (nn in 1:Nobs){
                yy <- nsucc[nn]
                mun <- mu[nn]
                dbinom1 <- dbinom(0:ntrials[nn], ntrials[nn],
                                  prob = mun)
                cumdbinom1 <- c(0,cumsum(dbinom1))
                lower[nn] <- cumdbinom1[yy + 1]
                upper[nn] <- cumdbinom1[yy + 2]
            }
        }
        if ((save.EQres == TRUE)|(start.control$selfstarts > 0)){
            unifEQres <- (lower + upper)/2
            ## Use inverse normal Phi^(-1)() to obtain normal EQres's:
            normEQres <- qnorm(unifEQres)
            ## Save to output:
            model.out$unifEQres <- unifEQres
            model.out$normEQres <- normEQres
        }
        if (save.RQres == TRUE){
            ## Unconditional RQres as the model has no clustering:
            unifs <- runif(Nobs)
            unifRQres <- lower + unifs * (upper - lower)
            ## Use inverse normal Phi^(-1)() to obtain normal RQRs:
            normRQres <- qnorm(unifRQres)
            ## Save to output:
            model.out$unifRQres <- unifRQres
            model.out$normRQres <- normRQres
        }
    }

    fixed.out <- model.out
    if (verbose > 1) cat("Finished unclustered part of model \n")
    
    ## Finished the glm() part of the model. 
    ## Stops here if no clustering. With clustering will overwrite model.out.
    ## Can use EQres and RQres values for starting clustering model fits.
        

    ## ---------------------------------------
    ## Do analysis if clustering is specified:
    ## ---------------------------------------

    ## Check specifications, get names sorted

    if (length(fact4clust) > 0){
        ## Set up some objects:
	Nlong <- Nobs*sum(nclust)
        pp.list <- NULL
        best.LLint <- -Inf
                 
        ## Number of factors for clustering
        nf4c <- length(fact4clust)
        if (nf4c > 2) stop("clustglm can cluster at most two factors")
        ## Set up the clustering:     
        if (verbose > 0) {
            cat("Setting up clustering of factor(s): \n")
            cat(fact4clust, "\n")
        }
        ## Error messages if not correctly specified
        if (length(nclust) != nf4c)
            stop("Please specify the number of clusters for each clustered factor")
        if (any(floor(nclust)==nclust)<1)
            stop("Please use integers for numbers of clusters")
        for (cf in 1:nf4c){
            if (!(fact4clust[cf] %in% names(data)))
                stop(paste("object",fact4clust[cf],"not found"))
            if (!is.factor(data[,names(data)==fact4clust[cf]]))
                stop(paste(fact4clust[cf],"is not a factor"))
        }

        ## Ask for starts if none given:
	## -----------------------------

	if (((start.control$selfstarts == 0)&(start.control$randstarts == 0))&
	    (is.null(start.control$alloc)))
            stop("Please choose an option in start.control")

	## Define objects used in the analysis:
	## ------------------------------------

        nlev <- rep(NA,nf4c)    ## Number of levels of factors to be clustered
        pp.list <- vector("list",nf4c)
        pi.list <- vector("list",nf4c)
        f4c.colno <- rep(NA,nf4c) ## Index of position of fact4clust in data frame
        for (cf in 1:nf4c) for (datacol in 1:ncol(data)){
            if (fact4clust[cf] == names(data)[datacol]){
                f4c.colno[cf] <- datacol
                nlev[cf] <- length(levels(data[,datacol]))
            }
        }
        if (verbose > 0){
            cat("Number of factors to be clustered = nf4c =", nf4c, "\n")
            cat("Number(s) of  clusters:", nclust, "\n")
            cat("Column number(s) in data frame = f4c.colno", "\n")
            cat(f4c.colno, "\n")
            cat("Number(s) of levels of factor(s) to be clustered = nlev", "\n")
            cat(nlev,"\n")
            cat("New name(s) of clustered factor(s):", "\n")
            cat(clustfactnames, "\n")
        }

        
        ## Stop if trivial clusterings:
        ## ----------------------------
        
        for (cn in 1:nf4c){          ## clustering number
            if (nclust[cn] == 1){
                cat("Trivial clustering of", fact4clust[cn], "into",
                    clustfactnames[cn], "clusters \n")
                cat("Only one cluster", "\n")
                stop("Please rewrite clustglm() inputs with only non-trivial clusterings")
            }
            if (nclust[cn] == nlev[cn]){
                cat("Trivial clustering with",fact4clust[cn],"the same as",
                    clustfactnames[cn],"\n")
                cat("Please replace",clustfactnames[cn],"with",
                    fact4clust[cn],"\n")
                stop("Please rewrite clustglm() inputs with only non-trivial clusterings")
            }
            if (nclust[cn] > nlev[cn]){
                stop("Number of clusters must be less than the number of levels of the factor")
            }
        }

        ## Will now have only non-trivial clusterings.

        
	## Check starting allocation if given:
        ## -----------------------------------
        
        ## Do error printout if alloc.start is given but not of right form
        ##    (checking: list, length of list, size of matrices)

        if (!is.null(start.control$alloc)){
	    alloc.start <- start.control$alloc
            if (verbose > 1) cat("Checking given starting allocation\n")
            if (is.matrix(alloc.start)){
                if (nf4c == 1) alloc.start <- list(pp.mat = alloc.start)
                if (nf4c == 2) 
                    stop("For biclustering, please supply two starting allocation matrices")
            }
            if (!is.list(alloc.start))
                stop("Please supply initial allocations to clusters as a list")
            if ((is.list(alloc.start))&(length(alloc.start) != nf4c))
                stop(paste("List of initial allocations should be of length",nf4c))
            if ((is.list(alloc.start))&(length(alloc.start) == nf4c)){
                for (cf in 1:nf4c){
                    if (!is.matrix(alloc.start[[cf]]))
                        stop(paste("Initial allocation of",fact4clust[cf],"to",
                                   clustfactnames[cf],"should be a matrix"))
                    if (is.matrix(alloc.start[[cf]])){
                        if (nrow(alloc.start[[cf]]) != nlev[cf])
                            stop(paste("Matrix of initial allocation of",fact4clust[cf],
                                       "to",clustfactnames[cf],"should have",nlev[cf],
                                       "rows"))
                        if (ncol(alloc.start[[cf]]) != nclust[cf])
                            stop(paste("Matrix of initial allocation of",fact4clust[cf],
                                       "to",clustfactnames[cf],"should have",
                                       nclust[cf],"columns"))
                    }
                }
            }
        }

        
	## -------------------------------
	## Set up selfstarts if requested:
	## -------------------------------

	if (start.control$selfstarts > 0){
	    if (verbose >0) cat("Creating SS.df","\n")
	    temp.df <- data
	    if (family == "gaussian") temp.df$Qres <- Qres
	    if ((family == "binomial")|(family == "poisson"))
                temp.df$Qres <- normEQres
	    ## Neg bin still to do
	    ## Use labels facA, facB, for factors in selfstart:
	    fo <- formula
	    inter.term <- attr(terms(fo),"term.labels")[(1:5)[attr(terms(fo),"order") == 2]]
	    strsp <- strsplit(inter.term,":")[[1]]
	    ## Might crash if formula has "*" not ":"
	    ## Look at first factor:
	    name1 <- strsp[1]
	    temp.df$facA <- temp.df[,names(temp.df) == name1]
	    nAclust <- length(levels(temp.df$facA))  ## Default if not clustering facA
	    if (strsp[1] %in% clustfactnames){
		name1 <- fact4clust[clustfactnames == strsp[1]]
		temp.df$facA <- temp.df[,names(temp.df) == name1]
		nAclust <- nclust[which(clustfactnames == strsp[1])]
            }
	    ## Look at second factor:
	    name2 <- strsp[2]
	    temp.df$facB <- temp.df[,names(temp.df) == name2]
	    nBclust <- length(levels(temp.df$facB))  ## Default if not clustering facB
	    if (strsp[2] %in% clustfactnames){
		name2 <- fact4clust[clustfactnames == strsp[2]]
		temp.df$facB <- temp.df[,names(temp.df) == name2]
		nBclust <- nclust[which(clustfactnames == strsp[2])]
	    }
	    ## Build SS.df for self-start.
	    SS.df <- data.frame(Qres = temp.df$Qres,
	                        facA = temp.df$facA, facB = temp.df$facB)
	    cat("Built SS.df")
        }

        
        ## --------------------------
        ## Construct long data frame:
        ## --------------------------

        if (verbose > 0) cat("Creating long data frame\n")
        
        ## For each clustering request, widen the data frame to include
        ## the new, clustered variables and lengthen data frame to allow
        ## for all combinations of levels of the new clustered variables.
        currentdata <- data
	currentdata$y.index <- as.factor(1:Nobs)

        cf.colno <- rep(NA,nf4c)      ## Column numbers of clustered factors.
        for (cf in 1:nf4c){
            ncl <- nclust[cf]
            longer.df <- currentdata
            if ((ncl > 1)&(ncl < nlev[cf]))
                for (repl in 2:ncl)
                    longer.df <- rbind(longer.df,currentdata)
            longer.df$x <- gl(n = ncl, k = nrow(currentdata),
                              length = nrow(longer.df),
                              labels = paste(clustfactnames[cf],seq_len(ncl),sep=""))
            names(longer.df)[ncol(longer.df)] <- clustfactnames[cf]
            cf.colno[cf] <- ncol(longer.df)
            currentdata <- longer.df
        }
        long.df <- longer.df

        rownames(long.df) <- as.character(1:nrow(long.df))

        if (verbose > 1) {
            cat("Structure of long.df","\n")
            cat(str(long.df),"\n")
            cat("Column number(s) of factor(s) to be clustered","\n")
            cat(f4c.colno,"\n")
            cat("Column number(s) for clustered factor(s)","\n")
            cat(cf.colno,"\n")
        }


        ## Find numbers of levels of all factors in long.df:
        ## -------------------------------------------------
        
        nlevels <- rep(NA,ncol(long.df))
        for (j in 1:ncol(long.df))
            if (is.factor(long.df[,j]))
                nlevels[j] <- length(levels(long.df[,j]))
        names(nlevels) <- names(long.df)

        if (verbose > 0) {
            cat("Numbers of levels of factors","\n")
            cat(nlevels,"\n")
        }
        
        
        ## Construct indicator matrices for the factors to be clustered:
        ## -------------------------------------------------------------
        
        ## ind[n,a] = 1 if row n has A at level a, else 0, n = 1 to N
        ## ind[n,b] = 1 if row n has B at level b, else 0, n = 1 to N
        
        indf4c.list <- vector("list",nf4c)
        for (f4c in 1:nf4c){
            indf4c.list[[f4c]] <- model.matrix(rpois(nrow(long.df),1) ~ -1 +
                                                   long.df[,f4c.colno[f4c]])
            dimnames(indf4c.list[[f4c]]) <- list(rownames(long.df),
                                                 levels(long.df[,f4c.colno[f4c]]))
        }
        names(indf4c.list) <- paste("indi",fact4clust,sep="")

        if (verbose > 0) cat("Indicator matrices done\n")

        if (verbose > 1) {
            cat("\n")
	    cat("Indicator matrix dimensions for factor(s) to be clustered\n")
	    cat(dim(indf4c.list[[1]]),"\n")
	    if (nf4c > 1) cat(dim(indf4c.list[[2]]),"\n")
        }


        ## Construct indicator matrices for the clustered factors:
        ## -------------------------------------------------------
        
        ## ind[n,u] = 1 if row n has U at level u, else 0, n = 1 to N
        ## ind[n,v] = 1 if row n has V at level v, else 0, n = 1 to N
        
        if (verbose > 1) cat("Column number(s) of clustered factor(s)",cf.colno,"\n")
        
        indcf.list <- vector("list",nf4c)
        for (cf in 1:nf4c){
            indcf.list[[cf]] <- model.matrix(rpois(nrow(long.df),1) ~ -1 +
                                                   long.df[,cf.colno[cf]])
            dimnames(indcf.list[[cf]]) <- list(rownames(long.df),
                                                 levels(long.df[,cf.colno[cf]]))
        }
        names(indcf.list) <- paste("indi",clustfactnames,sep="")
	
        if (verbose > 1) {
            cat("\n")
            cat("Indicator matrix dimensions for clustered factor(s)\n")
            cat(dim(indcf.list[[1]]),"\n")
            if (nf4c > 1) cat(dim(indcf.list[[2]]),"\n")
        }

	best.LLint <- -10^10
	    

        ## -----------------------
        ## Option 1: Random starts
        ## -----------------------
        
        ## If random starts chosen, for each do one M step,
	## choose the M-step output with the best logLik to initialise EM.

        randstarts <- start.control$randstarts
        if (randstarts>0){
            if (verbose > 0){
                cat("\n")
                cat("Trying", randstarts, "random starts\n")
            }

            for (rs in 1:randstarts){
                if (verbose > 0) cat("Random start number", rs, "\n")
                
                ## Do non-trivial random allocation:
                pp.list <- vector("list",nf4c)
                for (cf in 1:nf4c){
                    n1 <- nlev[cf]
                    n2 <- nclust[cf]
                    pp.list[[cf]] <- randpart(n1,n2)
                }
                if (rs == 1) best.pp <- pp.list

                ## Run one M step:
                
                M.out <- M_step(family, formula, response,
                                nf4c, nlev, nclust,
                                f4c.colno, cf.colno,
                                indcf.list, pp.list,
                                long.df, verbose)
                
                LLint <- M.out$LLint

		## If better than current best, save result.
                if (LLint > best.LLint){
		    best.M.out <- M.out
		    best.LLint <- LLint
		    if (verbose > 1) {
                        LLint.3 <- round(LLint,3)
                        cat("Found a better start. Best LLint so far:",LLint.3,"\n")
                    }
		}
		    
                ## Signal end of this random start:
                if (verbose > 0){
                    cat("Finished random start number ", rs, "\n")
                }	            
            
                ## Increment the loop:
                rs <- rs + 1
            }

            if (verbose > 0) cat("Finished", randstarts, " random starts\n")
        }

        ## Random starts done, best.M.out saved



        ## --------------------
        ## Option 2: Selfstarts
        ## --------------------

	SSdone <- "FALSE"
	
        if (start.control$selfstarts > 0){
	    no.selfstarts <- start.control$selfstarts
	    
            if (verbose > 0) {
                cat("\n")
                cat("Finding k-means using quantile residuals:\n")
            }
            ## Call the function for weighted single or double k-means:
	    ## The weighting is needed for unbalanced data

	    SS.alloc <- wsdkm(in.df = SS.df, nAclust = nAclust,
	                      nBclust = nBclust, randstart = no.selfstarts,
			      verbose)
            for (ii in 1:length(pp.list))
	        if (!all(dim(pp.list[[ii]]) == dim(SS.alloc[[ii]])))
		    stop("dimensions of pp matrices don't match")
            pp.list <- SS.alloc

            ## Do one M step using this pp.list:
            if (verbose > 1) cat("Starting the M step\n")
            this.glm <- fixed.glm
            M.out <- M_step(family, formula, response,
                            nf4c, nlev, nclust,
                            f4c.colno, cf.colno,
                            indcf.list, pp.list,
                            long.df, verbose)
            best.M.out <- M.out
            
            SSdone <- "TRUE"

            if (verbose > 0) cat("Finished selfstarts\n")
        }


        ## --------------------------------------------
        ## Option 3: Use a given starting allocation:
        ## --------------------------------------------
        
        ## Starting values for posterior probabilities:
        if (!is.null(start.control$alloc)){
	    if (SSdone == "TRUE")
	        warning("Given allocation not used as self start was successful")
	    alloc.start <- start.control$alloc
	    ## If single-mode clustering:
	    if (is.matrix(alloc.start)) pp.list <- list(alloc.start)
	    ## If biclustering:
            if (is.list(alloc.start)) pp.list <- alloc.start
            if (verbose > 1) {
                cat("\n")
                cat("Using the given starting allocation:\n")
            }
            if (verbose > 1) {
                pp1 <- round(pp.list[[1]],3)
                print(pp1)
                if (length(pp.list)==2) {
                    pp2 <- round(pp.list[[2]],3)
                    print(pp2)
                }
            }
            ## Do one M step:
            if (verbose > 1) cat("Starting an M step\n")
            this.glm <- fixed.glm
            M.out <- M_step(family, formula, response,
                            nf4c, nlev, nclust,
                            f4c.colno, cf.colno,
                            indcf.list, pp.list,
                            long.df, verbose)

            best.M.out <- M.out

            if (verbose > 0) cat("Finished start from given allocation\n")
        }

        ## -----------------------------------------------------------       
        ##                  Run the EM algorithm
        ## -----------------------------------------------------------       

        ## Initialise EM:
	
        options(warn=0)
        if (verbose > 0){
            cat("Setting up EM algorithm, doing up to ",
                EM.control$maxEMcycles," cycles","\n")
            cat("Starting from best M-step fit found so far\n")
        }

	M.out <- best.M.out
	LLint <- M.out$LLint
        LLc <- M.out$LLc
        long.df <- M.out$long.df
        this.glm <- M.out$this.glm
        pi.list <- M.out$pi.list
	## Need best.pp.list???
        
        ## Have just done an M step.
        ## Move to EM algorithm
        
        ## Make a vector of LL and parameters:
        EMoutvect <- c(LLint, LLc, this.glm$coef, c(pi.list[[1]]))
        if (length(pi.list) > 1) for (j in 2:length(pi.list))
            EMoutvect <- c(EMoutvect, c(pi.list[[j]]))
        ## Missing coefficients?
        EMoutvect[is.na(EMoutvect)] <- -10

        EMinvect  <- rep(1,length(EMoutvect))
        ## Save ests for printing?
        if (save.ests == TRUE){
            ests.seq <- matrix(NA, EM.control$maxEMcycles+1, length(EMoutvect)+1)
        }
        ## Do iterations:
        EMiter <- 1
      
        ## Run the EM cycles until the parameters have stabilised

        while(((EMiter==1)|(any(abs(EMinvect-EMoutvect)>EM.control$EMstoppingpar)))&
              (EMiter<(EM.control$maxEMcycles+1))){

            ## Do E step:
            pp.list <- E_step(nf4c, nlev, nclust,
                              f4c.colno, cf.colno,
                              pp.list, pi.list,
                              long.df, verbose)
            if (verbose >0){
                cat("EM iteration",EMiter,"\n")
                cat("E step done\n")
            }
            
            if (verbose > 1){
	        print(round(pp.list[[1]],3))
	        if (nf4c == 2) print(round(pp.list[[2]],3))
	    }

            ## Do M step:
            M.out <- M_step(family, formula, response,
                               nf4c, nlev, nclust,
                               f4c.colno, cf.colno,
                               indcf.list, pp.list,
                               long.df, verbose)
            pi.list <- M.out$pi.list    ## Proportions
            LLint <- M.out$LLint        ## Log likelihoods
            LLc <- M.out$LLc
            long.df <- M.out$long.df
            this.glm <- M.out$this.glm
            if (verbose > 1){
	        print(round(pi.list[[1]],3))
	        if (nf4c == 2) print(round(pi.list[[2]],3))
	    }

            if (verbose > 0){
                LLc.3 <- round(LLc,3)
                LLint.3 <- round(LLint,3)
                cat("M step done, LLc =", LLc.3, "LLint =", LLint.3,"\n")
            }

            ## Update vector of LL and pars, iteration number:
            EMinvect <- EMoutvect
            EMoutvect <- c(LLint, LLc, this.glm$coef,c(pi.list[[1]]))
            if (length(pi.list) == 2)
                EMoutvect <- c(EMoutvect,c(pi.list[[2]]))
            EMoutvect[is.na(EMoutvect)] <- -10
            if (save.ests == TRUE)
                ests.seq[EMiter,] <- c(EMiter, EMoutvect)
            EMiter <- EMiter + 1
        }        
        ## End of EM cycle

        if (EMiter<EM.control$maxEMcycles)
            if (verbose > 0) cat("EM algorithm converged after", EMiter-1, "iterations\n")
        if (EMiter>=(EM.control$maxEMcycles-1))
            if (verbose > 0) cat("EM algorithm stopped after", EMiter-1, "iterations\n")

        
        ## Save results for output:
	## ------------------------
	
        this.glm <- M.out$this.glm
        LLc <- M.out$LLc
        LLint <- M.out$LLint

        if (save.ests == TRUE){
            colnames.ests <- c("iter", "LLint", "LLc",
                               paste("coef", 1:length(this.glm$coef), sep = ""),
                               paste("pi", 1:length(pi.list[[1]]), sep = ""))
            if (length(pi.list) == 2) colnames.ests <-
                c(colnames.ests, paste("kappa", 1:length(pi.list[[2]]), sep = ""))
            colnames(ests.seq) <- colnames.ests
        }

        ## Warning messages
        ## Failure to distinguish enough groups:
        if (verbose > 0)  {
            for (j in 1:length(pi.list)) {
                ## pi.list[[j]][is.na(pi.list[[j]])] <- 0
                if (min(pi.list[[j]], na.rm=T)==0)
                    cat("Cannot distinguish", nclust[j],
                        "groups for Factor", clustfactnames[j],
                        ". Try re-running the model.\n")
            }
        }

        ## Collect terms which match the unclustered glm output:

        model.out <- list(coefficients = this.glm$coefficients,
                          residuals = this.glm$residuals,
                          fitted.values = this.glm$fitted.values,
                          rank = this.glm$rank,
                          family = this.glm$family,
                          linear.predictors = this.glm$linear.predictors,
                          deviance = this.glm$deviance,
                          null.deviance = this.glm$null.deviance,
                          weights = this.glm$weights,
                          df.residuals = this.glm$df.residual,
                          df.null = this.glm$df.null,
                          y = this.glm$y,
                          model = this.glm$model,
                          call = this.glm$call,
                          formula = this.glm$formula,
                          terms = this.glm$terms,
                          data.in = data,
                          data = this.glm$data)
 
        ## Overwrite terms in model.out for a clustglm object:
        model.out$formula <- formula
        model.out$deviance <- -2*LLint

        ## Add terms linked to clustering:
        npar <- model.out$rank + sum(nclust) - length(nclust)
        model.out$npar <- npar
        model.out$LLint <- LLint
        model.out$LLc <- LLc
        model.out$AIC <- -2*LLint + 2*npar
        model.out$BIC <- -2*LLint + npar*log(Nobs)
        model.out$fact4clust <- fact4clust
        model.out$clustfactnames <- clustfactnames
        model.out$nclust <- nclust
        model.out$data.in <- data
        model.out$start.control <- start.control
        model.out$EM.control <- EM.control
        model.out$EMcycles <- EMiter - 1
	names(pi.list) <- clustfactnames
        for (cf in 1:nf4c) {
            names(pi.list[[cf]]) <- levels(long.df[,cf.colno[cf]])
        }
        model.out$pi.list <- pi.list
        names(pp.list) <- clustfactnames
        for (cf in 1:nf4c) {
            dimnames(pp.list[[cf]]) <- list(levels(long.df[,f4c.colno[cf]]), 
                                            levels(long.df[,cf.colno[cf]]))
        }
        model.out$pp.list <- pp.list
        model.out$final.glm <- this.glm
	
        ## Add unconditional fitted values (ufits) to output:
	yi <- long.df$y.index
	ufits <- rep(NA,Nobs)
	for (ii in 1:Nobs){
	    tiny.df <- long.df[yi == levels(yi)[ii],]
	    ufits[ii] <- sum(tiny.df$wts * tiny.df$mu)
	}
	model.out$uncond.fits <- ufits

        ## Ancillary parameters:
	
        if (family[1] == "gaussian")
            model.out$sigsq <- long.df$sigsq[1]

        ## Optional saved items:
        
        if (save.long == TRUE){
            ## names(long.df)[names(long.df) == "mu"] <- "cond.fit"
            model.out$long.df <- long.df
        }
    
        if (save.ests == TRUE){
	    ## Remove the first row of the matrix, used in setting up:
            ests.seq <- ests.seq[!is.na(ests.seq[,1]),]
        model.out$ests.seq <- ests.seq
        }
	
        if (save.Qres == TRUE){
            ## Qres for continuous data is a PIT (prob. integral transform).
            ## Only implemented for Gaussian so far.
            if (family == "gaussian"){
		## Set up conditional Qres vector, one per row of long dataframe:
                condQres <- rep(NA,Nlong)
		## FOLLOWING NEEDS CHECKING
                mu <- this.glm$fitted.values  ## Conditional fitted values
                y <- this.glm$y
                sigma <- sqrt(long.df$sigsq[1])  ## Check, sigsq from M step
                for (nn in 1:Nlong){
                    yy <- y[nn]
                    mun <- mu[nn]
                    condQres[nn] <- pnorm(yy, mean = mun, sd = sigma)
                }
		long.df$condQres <- condQres
                ## Other continuous distributions still to be implemented
                ## Need one quantile residual per data point.
                ## Make the conditional QR unconditional
                this.df <- this.glm$data
                yi <- this.df$y.index
                unifQres <- rep(NA,Nobs)
                for (ii in 1:Nobs){
                    tiny.df <- long.df[yi == levels(yi)[ii],]
                    unifQres[ii] <- sum(tiny.df$wts * tiny.df$condQres)
                }
		model.out$unifQres <- unifQres
                ## Use inverse normal Phi^(-1)() to obtain normal Qres values:
                normQres <- qnorm(unifQres)
		model.out$normQres <- normQres
            }
        }
  
        ## Now do discrete data, need lower and upper ends of risers.
        if ((save.RQres == TRUE)|(save.EQres == TRUE)){
            if ((family == "binomial")|(family == "poisson")){
                ## Set up lower and upper limits of risers:
                lower <- rep(NA,Nlong)
                upper <- rep(NA,Nlong)
                mu <- long.df$mu
                ## Assuming Poisson family
                if (family == "poisson"){
                    y <- long.df$y
                    MM <- max(y)
                    for (nn in 1:Nlong){
                        yy <- y[nn]
                        mun <- mu[nn]
                        dpois1 <- dpois(0:MM, lambda = mun)
                        cumdpois1 <- c(0,cumsum(dpois1),1)
                        lower[nn] <- cumdpois1[yy+1]
                        upper[nn] <- cumdpois1[yy+2]
                    }
                }
                ## Assuming binomial family:
                if (family == "binomial"){
                    nsucc <- model.out$data$nsucc
                    ntrials <- model.out$data$ntrials
                    for (nn in 1:Nlong){
                        yy <- nsucc[nn]
                        mun <- mu[nn]
                        dbinom1 <- dbinom(0:ntrials[nn], ntrials[nn],
                                          prob = mun)
                        cumdbinom1 <- c(0,cumsum(dbinom1))
                        lower[nn] <- cumdbinom1[yy + 1]
                        upper[nn] <- cumdbinom1[yy + 2]
                    }
                }
		## Find conditional and unconditional RQres values:
		if (save.RQres == TRUE){
                    ## For each original data point, assign one runif[0,1] value,
		    ## save in long vector.
		    runifs <- rep(runif(Nobs), sum(nclust))
		    ## Find conditional RQres values, save in long.df.
		    long.df$condRQres <- lower + runifs * (upper - lower)
		    ## Find unconditional RQres values, one per datum.
		    ## Save both uniform and normal versions.
		    yi <- long.df$y.index
		    unifRQres <- rep(NA,Nobs)
		    for (n in 1:Nobs){
		        tiny.df <- long.df[yi == levels(yi)[n],]
                        unifRQres[n] <- (sum(tiny.df$wts * tiny.df$condRQres))/
			                (sum(tiny.df$wts))
                    }
		    ## Use inverse normal Phi^(-1)() to obtain normal RQres:
                    normRQres <- qnorm(unifRQres)
                    ## Save to output:
                    model.out$unifRQres <- unifRQres
                    model.out$normRQres <- normRQres
               }

	       ## Find conditional and unconditional EQres values:
	       if (save.EQres == TRUE){
	           long.df$condEQres <- (upper + lower)/2   ## Length Nlong                
		    ## Find unconditional EQres values, one per datum.
		    ## Save both uniform and normal versions.
		    yi <- long.df$y.index
		    unifEQres <- rep(NA,Nobs)
		    for (n in 1:Nobs){
		        tiny.df <- long.df[yi == levels(yi)[n],]
                        unifEQres[n] <- (sum(tiny.df$wts * tiny.df$condEQres))/
			                (sum(tiny.df$wts))
                    }
		    ## Use inverse normal Phi^(-1)() to obtain normal EQres:
                    normEQres <- qnorm(unifEQres)
                    ## Save to output:
                    model.out$unifEQres <- unifEQres
                    model.out$normEQres <- normEQres
               }
           }
        }
    ## Save fixed-part model:
    model.out$fixed.part.model <- fixed.out
    }
    ## Finished model.out for clustered model.
    
    class(model.out) <- c('clustglm','glm')
    
    return(model.out)
}
