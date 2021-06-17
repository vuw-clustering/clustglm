## clustbiplot.R


## The clustglm version of a biplot.
## Assume biclust 3 by 3 model, data from a matrix, one obs. per cell.
## Assume at most two main effects terms, row and col from data matrix,
## to be labelled internally as factors Af, Bf resp.
## Assume one interaction term, both clustered factors, labelled
## internally as Uf and Vf resp.
## Need three gamma matrices, gammaUV for biclustering, gammaUB for
## row clustering and gammaAV for col clustering.

clustbiplot <- function(model,
                        data.mat,
                        verbose = 1,
                        ## Graphical parameters:
			main = NULL, sub = NULL,
                        xlim = NULL, ylim = NULL,
                        xlab = "Dimension 1", ylab = "Dimension 2",
                        pch = NULL,
                        rnames = NULL, cnames = NULL,
                        rcnames = NULL, ccnames = NULL,
                        colrvect = NULL,
                        myratio = 1,
                        centroids = TRUE,
                        partitions = TRUE,
                        axis.lines = TRUE){
    
    this.model <- model
## this.model <- bi33.cglm
## data.mat <- counts.mat

    partitions = TRUE
    axis.lines = TRUE

    ## Y matrix:
    Y.m <- data.mat
    
    ## Estimated Z and X matrices:   
    Z.m <- this.model$pp.list[[1]]       ## A by U
    X.m <- this.model$pp.list[[2]]       ## B by V

    ## Estimated pi and kappa vectors:
    pi.v <- this.model$pi.list[[1]]
    kappa.v <- this.model$pi.list[[2]]

    ## Give chosen names to levels of factors:
    if (is.null(rnames)) rnames <- rownames(Y.m)
    if (is.null(cnames)) cnames <- colnames(Y.m)

    ## Use Z and X matrices to do psi estimates (PA2014 eq. 11):
    ## ---------------------------------------------------------

    Tot <- sum(Y.m)
    Ta <- apply(Y.m,1,sum)
    Tb <- apply(Y.m,2,sum)
    ZTYX <- t(Z.m) %*% Y.m %*% X.m  ## Matrix, U=3 by V=3
    YX <- Y.m %*% X.m     ## Matrix, A by 3 
    ZTY <- t(Z.m) %*% Y.m     ## Matrix, 3 by B 
    ZTya <- t(Z.m) %*% Ta     ## Vector, length U=3
    ybX <- Tb %*% X.m         ## Vector, length V=3

    psiUV <- (Tot * ZTYX) / (ZTya)%*%(ybX) 
    deltaUV <- psiUV - 1

    psiUB <- (Tot * ZTY) / (ZTya)%*%(Tb) 
    deltaUB <- psiUB - 1

    psiAV <- (Tot * YX) / (Ta)%*%(ybX) 
    deltaAV <- psiAV - 1



    ## Do the plot:
    ## ------------

    if (is.character(colrvect)){
        if (length(colrvect) == 1){
            colrvect <- rep("colrvect",6)
	}
    }
    if (is.null(colrvect))
        colrvect <- c("red","blue","purple","green","orange2","tan3")

    ## Do matrix calculations, SVD of gamma:
    svdUV <- svd(deltaUV)
    Usvd <- svdUV$u
    Vsvd <- svdUV$v
    d.vect <- svdUV$d  ## To be partitioned equally using sqrt
    d.vect[d.vect<0] <- 0  ## In case of rounding problems
    dnames <- paste("D",1:length(d.vect),sep="")  ## Dimension names
    dimnames(Usvd) <- list(rownames(deltaUV),dnames)
    dimnames(Vsvd) <- list(colnames(deltaUV),dnames)
    S.mat <- diag(d.vect)
    dimnames(S.mat) <- list(dnames,dnames)
    myrat <- c(myratio/2,1-myratio/2)
    S1.mat <- diag(d.vect^myrat[1])
    S2.mat <- diag(d.vect^myrat[2])
    dimnames(S1.mat) <- list(dnames,dnames)
    dimnames(S2.mat) <- list(dnames,dnames)
    
    ## Find points for cluster centroids:
    G.mat <- Usvd %*% S1.mat     ## Matrix G (U by 3)
    rowclustpts <- G.mat
    H.mat <- Vsvd %*% S2.mat     ## Matrix H (V by 3)
    colclustpts <- H.mat

    ## Find associated row points and column points:
    ## Use same operations as for centroids
    S1.mat[3,3] <- max(S1.mat[3,3],1e-20)
    S2.mat[3,3] <- max(S2.mat[3,3],1e-20)
    rowpts <- deltaAV %*% Vsvd %*% solve(S1.mat)
    colpts <- t(deltaUB) %*% Usvd %*% solve(S2.mat)

    ## Create the plot:
    all.mat <- rbind(G.mat,H.mat,rowpts,colpts)

    plot(all.mat[,1:2],axes=F,type="n",
             xlab = xlab, ylab = ylab,las=1,
             main = main)
    axis(1)
    axis(2)
    box()
    rnames <- rownames(rowpts)
    cnames <- rownames(colpts)
    text(rowpts[,1:2], rnames, col=colrvect[1])
    text(colpts[,1:2], cnames, col=colrvect[2])
    if (centroids == TRUE){
        rclustnames <- paste("rc",1:3,sep="")
        cclustnames <- paste("cc",1:3,sep="")
        text(rowclustpts[,1:2], rclustnames, cex = 1.5, col=colrvect[1])
        text(colclustpts[,1:2], cclustnames, cex = 1.5, col=colrvect[2])
    }


    ## pdf(file = "count.biplot.pdf")

    plot(all.mat[,1:2],axes=F,type="n",
             xlab = xlab, ylab = ylab,las=1,
             main = main)
    axis(1)
    axis(2)
    box()
    rnames <- rownames(rowpts)
    cnames <- rownames(colpts)
    text(rowpts[,1:2], rnames, col=colrvect[1])
    text(colpts[,1:2], cnames, col=colrvect[2])
    
    if (centroids == TRUE){
        rclustnames <- paste("rc",1:3,sep="")
        cclustnames <- paste("cc",1:3,sep="")
        points(rowclustpts[,1:2], pch=10, cex = 1.5, col=colrvect[1])
        points(colclustpts[,1:2], pch=10, cex = 1.5, col=colrvect[2])
        for (i in 1:3){
            lines(c(0,rowclustpts[i,1]),
                  c(0,rowclustpts[i,2]),
                  col=colrvect[1])
        }
        for (i in 1:3){
            lines(c(0,colclustpts[i,1]),
                  c(0,colclustpts[i,2]),
                  col=colrvect[2])
        }
    }

    if (axis.lines == TRUE){
        lines(c(2*min(all.mat[,1]),2*max(all.mat[,1])),c(0,0),
              col="grey",lty=5)
        lines(c(0,0),c(2*min(all.mat[,2]),2*max(all.mat[,2])),
              col="grey",lty=5)
    }

    if (partitions == TRUE){
        cat("partitions not yet implemented\n")
    }

    ## dev.off()

    ## Return all.mat as a df:
    all.df <- as.data.frame(all.mat)
    return(all.df)
}




if (1 == 0){

    library(MASS)

## pdf(file = "count.biplot.eqsc.pdf")

    eqscplot(all.mat[,1:2],axes=F,type="n",
             xlab = xlab, ylab = ylab,las=1,
             main = main, ratio = myratio)
    axis(1)
    axis(2)
    box()
    rnames <- rownames(rowpts)
    cnames <- rownames(colpts)
    text(rowpts[,1:2], rnames, col=colrvect[1])
    text(colpts[,1:2], cnames, col=colrvect[2])
    if (centroids == TRUE){
        rclustnames <- paste("rc",1:3,sep="")
        cclustnames <- paste("cc",1:3,sep="")
        text(rowclustpts[,1:2], rclustnames, cex = 1.5, col=colrvect[1])
        text(colclustpts[,1:2], cclustnames, cex = 1.5, col=colrvect[2])
    }

## dev.off()
}
