## Function findstart()

## For matrix data, use other clustering methods to try to find a good start.
## Use if the model has the clustered factor(s), rowclust and/or colclust, occurring
## only in an interaction term.
## (i)   Fit the known, unclustered part of the model.
## (ii)  Use the expected quantile residuals from (i) with k-means or double k-means
##       to obtain a preliminary clustering of the rows and/or columns.
## (iii) Use the clusters in (ii) for a starting allocation in clustglm().

findstart <- function(x, clustmode = "rows", nclust,
                      clusttype = "kmeans", randstarts = 1,
                      metric = "euclidean", dendro = FALSE){

    if(!is.matrix(x)) stop("x needs to be a data matrix")

    if (clustmode == "cols") x <- t(x)

    if (!(clustmode == "both")){
        ## Single mode clustering
        
        if (clusttype == "kmeans"){
            km <- kmeans(x, centers = nclust[1], nstart = randstarts)
            pp1 <- matrix(0,length(km$cluster),nrow(km$centers))
            rownames(pp1) <- names(km$cluster)
            colnames(pp1) <- rownames(km$centers)
            for (clustno in 1:ncol(pp1))
                pp1[km$cluster == clustno,clustno] <- 1
            pp.start <- list(pp1)
        }

        if (clusttype == "hclust"){
            hc <- hclust(d = dist(x, method = metric))
            if (dendro == TRUE) plot(hc)
            hc.cut <- cutree(hc, k = nclust[1])
            pp1 <- matrix(0,length(hc.cut),nclust[1])
            rownames(pp1) <- names(hc.cut)
            colnames(pp1) <- as.character(1:nclust[1])
            for (clustno in 1:ncol(pp1))
                pp1[hc.cut == clustno,clustno] <- 1
            pp.start <- list(pp1)
        }

        if (clusttype == "diana"){
            di <- diana(x = dist(x, method = metric))
            if (dendro == TRUE) plot(di)
            di.cut <- cutree(di, k = nclust[1])
            pp1 <- matrix(0,length(di.cut),nclust[1])
            rownames(pp1) <- rownames(x)
            colnames(pp1) <- as.character(1:nclust[1])
            for (clustno in 1:ncol(pp1))
                pp1[di.cut == clustno,clustno] <- 1
            pp.start <- list(pp1)
        }

    }

    if (clustmode == "both"){
        if (clusttype != "dkm")
            warning("Using double kmeans (dkm), since this is biclustering")
        if (length(nclust) != 2)
            stop(" Need the input nclust to be of length two for biclustering")
        dkm.out <- dkm(y.matrix = x,
                       rowclust = nclust[1],
                       colclust = nclust[2],
                       randstart = randstarts)
        pp.start <- list(pp1 = dkm.out[[1]],pp2 = dkm.out[[2]])
        rownames(pp.start[[1]]) <- rownames(x)
        colnames(pp.start[[1]]) <- as.character(1:nclust[1])
        rownames(pp.start[[2]]) <- colnames(x)
        colnames(pp.start[[2]]) <- as.character(1:nclust[2])
    }

    pp.start
}



