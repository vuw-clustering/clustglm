## ordplot.R

## Construct ordpoints.mat, e.g. on hyperplane,
## for all c, sum pi(r) gamma(rc) = 0.
## Use SVD to reduce one dimension, e.g. plane in 3D rotated to
## have 3rd coordinate zero, so in 2D.


ordplot <- function(model, ord.factor, clust.factor,
                    main = NULL, sub = NULL,
                    xlim = NULL, ylim = NULL,
                    xlab = "SVD1", ylab = "SVD2",
                    pch = 4,
                    pt.labels = NULL, 
                    colourvector = NULL){
    this.model <- model
    long.df <- this.model$long.df
    fitted.v <- this.model$fitted
    
    ## Find values of factors:
    ord.f <- long.df[,(1:ncol(long.df))[ord.factor == names(long.df)]]
    clust.f <- long.df[,(1:ncol(long.df))[clust.factor == names(long.df)]]
    nx <- length(unique(ord.f))
    ncl <- length(unique(clust.f))
    if (this.model$family$family == "poisson")
        if (!(ncl == 3))
            stop("For Poisson family data use three clusters")
    if (!(this.model$family$family == "poisson"))
        if (!(ncl == 2))
            stop("For non-Poisson data use two clusters")
    
    ## Find the proportions vector and gamma matrix for plotting:
    ## findpars() will also check that there is one interaction term.
    these.pars <- findpars(this.model)
    gam.mat <- these.pars[[length(these.pars)]]
    ## Needs to be ord.factor by clust.factor
    if (all(rownames(gam.mat)[1:2] == levels(clust.f)[1:2]))
        gam.mat <- t(gam.mat)
    if (this.model$family$family == "poisson"){
        psi.mat <- exp(gam.mat)
        pp.mat <- this.model$pp.list[[1]]
        prop.vec <- apply(pp.mat,2,mean)
        ## print(prop.vec)
        ## Calculate matrix of ordination points:
        ordpoints.mat <- psi.mat %*% diag(prop.vec)
        ## print(ordpoints.mat)
        if (is.null(pt.labels)) pt.labels <- rownames(ordpoints.mat)
        ## print(pt.labels)
        ## Have three clusters
        ## Supply a main title:
        if (is.null(main)) main =
            paste("2D ordination of",ord.factor,"using",ncl,clust.factor)
        this.svd <- svd(ordpoints.mat)
        this.svd$d[this.svd$d < 0.000001] <- 0.000001
        plotpts.mat <- this.svd$u %*% diag(sqrt(this.svd$d))
        rownames(plotpts.mat) <- pt.labels
        ## print(plotpts.mat)
        plot(plotpts.mat[,1:2], type = "n",
             xlab = xlab, ylab = ylab,
             main = main)
        lines(range(plotpts.mat[,1]),c(0,0),lty=5,col="grey")
        lines(c(0,0),range(plotpts.mat[,2]),lty=5,col="grey")
        text(plotpts.mat[,1:2],pt.labels,cex=1.3)
    }
    if (!(this.model$family$family == "poisson")){
        ## Have ncl = 2
        stop("Ordplot not yet implemented for non-Poisson data")
        ## Supply a main title:
        ## if (is.null(main)) main =
        ##     paste("1D ordination of",ord.factor,"using",ncl,clust.factor)
        ## this.svd <- svd(ordpoints.mat)
        ## this.svd$d[this.svd$d < 0.000001] <- 0.000001
        ## plotpts.mat <- this.svd$u %*% diag(sqrt(this.svd$d))
        ## rownames(plotpts.mat) <- pt.labels
        ## ## print(plotpts.mat)
        ## plot(c(-1,1),range(plotpts.mat[,1]), type = 'n',
        ##      axes = F, xlab = "", ylab = xlab,
        ##      main = main)
        ## axis(2)
        ## box()
        ## lines(c(0,0),range(plotpts.mat[,1]))
        ## lines(c(-1,1),c(0,0),lty=5,col="grey")
        ## points(rep(0,nrow(plotpts.mat)),plotpts.mat[,1],
        ##        pch = pch)
        ## text(rep(0.1,nrow(plotpts.mat)),plotpts.mat[,1],
        ##      pt.labels,cex=1.3)
    }
    

    ## Return the data frame:
    pts.df <- as.data.frame(plotpts.mat)
    return(pts.df)
}

