## comparison.R

comparison <- function(modellist,rounding = 3){

    ## Construct a data frame which tables LL, ResDev, npar, AIC, BIC, relAIC, relBIC
    ## by model.
    M <-length(modellist)
    this.mat <- matrix(NA,M,7)
    rownames(this.mat) <- names(modellist)
    colnames(this.mat) <- c("logL","Res.Dev.","npar","AIC","relAIC","BIC","relBIC")
    for (m in 1:M){
        this.mod <- modellist[[m]]
	if (is.numeric(this.mod$LL))  ## Unclustered model
            this.mat[m,1] <- this.mod$LL
	if (is.numeric(this.mod$LLint))  ## Clustered model
            this.mat[m,1] <- this.mod$LLint
        this.mat[m,3] <- this.mod$npar
        this.mat[m,4] <- this.mod$AIC
        this.mat[m,6] <- this.mod$BIC
    }
    this.mat[,2] <- -2*this.mat[,1]
    this.mat[,5] <- this.mat[,4] - min(this.mat[,4])
    this.mat[,7] <- this.mat[,6] - min(this.mat[,6])
    ## Print to screen and save if required.
    print(round(this.mat,rounding))
    out.df <- as.data.frame(this.mat)
    invisible(out.df)
}
