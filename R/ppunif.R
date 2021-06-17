## ppunif.R

## Diagnostic PPplot using randomised quantile residuals, RQR.
## If model fits data well, a plot of RQR versus Unif{0,1}
## quantiles should be a straight line.
## This is a weighted PP plot, to allow for some "points"
## having weights less than one, because of the mixture
## models.


ppunif <- function(model, main = NULL, uline = TRUE,
                   xlab = NULL, ylab = NULL){
    this.model <- model
    if (is.null(xlab)) xlab <- "Uniform[0,1] cdf"
    if (is.null(ylab)) ylab <- "Quantile residual"
    ## Case 1, not a clustered model:
    if (is.null(this.model$pp.list)){
        RQR <- sort(this.model$RandQuantRes)
        wtRQR <- RQR
        n <- length(wtRQR)
        UQuants <- ((1:n)-0.5)/n
    }
    ## Case 2, clustered model:
    if (is.list(this.model$pp.list)){
        ## Start with conditional RQRs:
        RQR <- this.model$RandQuantRes
        wts <- this.model$data$wts
        long <- this.model$data
        wtRQR <- wts*RQR
        long$wtRQR <- wtRQR
        data.in <- this.model$data.in
        
        ## For each y create a single RQR,
        ## using total probability formula.
        oneRQR <- tapply(wtRQR,long$y.index,sum)
        
        ## Sort and construct plotting values:
        data.out <- data.in
        data.out$RQR <- oneRQR
        wtRQR <- sort(oneRQR)
        n <- length(wtRQR)
        UQuants <- ((1:n)-0.5)/n
    }

    ## Do ppunif plot:
    plot(UQuants, wtRQR, pch = 16, main=main,
         xlim = c(0,1), ylim = c(0,1),
         xlab = xlab, ylab = ylab)
    if (uline) lines(c(0,1),c(0,1))
    
    ## Return the data frame invisibly:
    ppunif.df <- data.frame(Uniform.quants  = UQuants,
                            WtedRQR = wtRQR)
    invisible(ppunif.df)
}
