## qqunif.R

qqunif <- function(y, main = "Uniform Q-Q Plot",
                   xlab = "Theoretical Quantiles",
                   ylab = "Sample Quantiles"){
    n <- length(y)
    unifquants <- ((1:n)-0.5)/n
    plot(unifquants, sort(y),
         main = main,
         xlab = xlab,
         ylab = ylab)
}

## qqunifline <- function(y, distribution = qunif,
##                       probs = c(0.25, 0.75),
## not implemented

    
