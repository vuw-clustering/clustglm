## findQR.R

## Find quantile residuals (QRs) for a model fitted to continuous data.
## Currently only implemented for gaussian data.

findQR <- function(model){
    this.model <- model
    ## Check it is a clustglm object:
    if (class(this.model)[1] != "clustglm")
        stop("The function findQR can only be used on a model of class 'clustglm'.")
    if (this.model$family[[1]] != "gaussian")
                cat("No quantile residual calculated, family not gaussian.\n")
    ## With continuous data, the QR is just a probability integral transform:
    if (this.model$family[1] == "gaussian"){
        ## QR is just a probability integral transform:
        NN <- nrow(this.model$data)
        condQR <- rep(NA,NN)
        mu <- this.model$fitted.values
        y <- this.model$y
        sigma <- sqrt(this.model$sigsq[1])
        for (nn in 1:NN){
            yy <- y[nn]
            mun <- mu[nn]
            condQR[nn] <- pnorm(yy, mean = mun, sd = sigma)
        }
        ## Other distributions still to be implemented
        
        ## Need one quantile residual per data point.
        ## If it is a clustered model:
        if (is.numeric(this.model$nclus)){
            ## Make the conditional QR unconditional
            this.df <- this.model$data
            yi <- this.df$y.index
            lyi <- length(levels(yi))
            unifQR <- rep(NA,lyi)
            for (ii in 1:lyi){
                tiny.df <- this.df[yi == levels(yi)[ii],]
                unifQR[ii] <- sum(tiny.df$wts * tiny.df$condQR)
            }
        }
        ## If the model has no clustering:
        if (is.null(this.model$nclus)){
            unifQR <- condQR
        }
                
        ## Now use inverse normal Phi^(-1)() to obtain normal QRs:
        normQR <- qnorm(unifQR)
        
        ## Return a data frame with unifQR and normQR:
        QR <- data.frame(unifQR = unifQR, normQR = normQR)
        if (is.numeric(this.model$nclus)){
            rownames(QR) <- levels(yi)
        }

    }
    return(QR)
}
