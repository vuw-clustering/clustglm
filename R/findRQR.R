## findRQR.R

## Generate randomised quantile residuals (RQRs) for a fitted model.
## For the same model, RQRs will differ every time.
## If the model fits the data well, the RQRs will follow a N(0,1)
## distribution. 

findRQR <- function(model){
    this.model <- model
    ## Check it is a clustglm object:
    if (class(this.model)[1] != "clustglm")
        stop("The function findRQR can only be used on a model of class 'clustglm'.")
    if ((this.model$family[[1]] != "poisson") & (this.model$family[[1]] != "binomial"))
        stop("Try findQR; the function findRQR can only be used on discrete data")
    ## Set up lower and upper limits of risers:
    NN <- nrow(this.model$data)
    lower <- rep(NA,NN)
    upper <- rep(NA,NN)
    mu <- this.model$fitted.values
    ## Assuming Poisson family
    if (this.model$family[[1]] == "poisson"){
        y <- this.model$y
        MM <- max(y)
        for (nn in 1:NN){
            yy <- y[nn]
            mun <- mu[nn]
            dpois1 <- dpois(0:MM, lambda = mun)
            cumdpois1 <- c(0,cumsum(dpois1),1)
            lower[nn] <- cumdpois1[yy+1]
            upper[nn] <- cumdpois1[yy+2]
        }
    }
    ## Assuming binomial family:
    if (this.model$family[[1]] == "binomial"){
        nsucc <- this.model$data$nsucc
        ntrials <- this.model$data$ntrials
        for (nn in 1:NN){
            yy <- nsucc[nn]
            mun <- mu[nn]
            dbinom1 <- dbinom(0:ntrials[nn], ntrials[nn],
                              prob = mun)
            cumdbinom1 <- c(0,cumsum(dbinom1))
            lower[nn] <- cumdbinom1[yy + 1]
            upper[nn] <- cumdbinom1[yy + 2]
        }
    }

    ## For each original data point, assign one runif[0,1] value:
    ## If it is a clustered model:
    if (is.numeric(this.model$nclus)){
        this.df <- this.model$data
        yi <- this.df$y.index
	lyi <- length(levels(yi))
	unifs <- rep(NA,NN)
	for (ii in 1:lyi){
 	   this.runif <- runif(1)
    	   unifs[yi == levels(yi)[ii]] <- this.runif
        }
        ## Find conditional unifRQRs:
        this.df$condUnifRQR <- lower + unifs * (upper- lower)
	## Each of the above conditional unifRQR values must be made
        ## unconditional, giving one value for each original datum.
        unifRQR <- rep(NA,lyi)
        for (ii in 1:lyi){
            tiny.df <- this.df[yi == levels(yi)[ii],]
    	    unifRQR[ii] <- sum(tiny.df$wts * tiny.df$condUnifRQR)
        }
    }
    ## If the model has no clustering:
    if (is.null(this.model$nclus)){
	unifs <- runif(NN)
	unifRQR <- lower + unifs * (upper - lower)
    }
    ## Now use inverse normal Phi^(-1)() to obtain normal RQRs:
    normRQR <- qnorm(unifRQR)

    ## Return a data frame with unifRQR and normRQR:
    RQR <- data.frame(unifRQR = unifRQR, normRQR = normRQR)
    if (is.numeric(this.model$nclus)){
        rownames(RQR) <- levels(yi)
    }
    
    return(RQR)
}


