## dkm.R

## ---------------
## Double k-Means:
## ---------------

dkm <- function(y.matrix, rowclust, colclust, randstart){
    
    maxiter <- 100
    eps <- 0.0000000001   ## convergence tolerance
    n <- nrow(y.matrix)
    p <- ncol(y.matrix)
    ys.mat <- y.matrix         ## Standardised data matrix
    st <- sum(ys.mat^2)   ## Overall SS:
    
    for (loop in 1:randstart){
        ## generate a random partition for variables (cols)
        ppc.ma <- randpart(p,colclust)
        ## generate a random partition for units     (rows) 
        ppr.ma <- randpart(n,rowclust)
        ## generate matrix of means
        su <- apply(ppr.ma,2,sum)  # No. rows in rowclust
        sv <- apply(ppc.ma,2,sum)  # No. cols in colclust
        Ym <- diag(1/su) %*% t(ppr.ma) %*% ys.mat %*% ppc.ma %*% diag(1/sv)
        ## Ym is rowclust by colclust, mean values in each rowclust and colclust
        B <- ppr.ma %*% Ym %*% t(ppc.ma)   # n by p, fitted values
        f0 <- sum(diag((t(B) %*% B)))/st   # explained variance
        it <- 0
        
        ## iteration phase, update the unknown allocations
        while (it<=maxiter){
            it <- it+1
            ## given ys.mat and ppc.ma update ppr.ma
            ppr.ma <- matrix(0,n,rowclust)
            Ymv <- Ym %*% t(ppc.ma)   # Ymv is rowclust by p
            for (i in 1:n){
                mindif <- sum((ys.mat[i,]-Ymv[1,])^2)
                posmin <- 1
                for (k in 2:rowclust){
                    dif <- sum((ys.mat[i,]-Ymv[k,])^2)
                    if (dif < mindif){
                        mindif <- dif
                        posmin <- k
                    } 
                }
                ppr.ma[i,posmin] <- 1
            }
            su <- apply(ppr.ma,2,sum)
            
            ## if a class is empty split the largest class
            while (min(su)==0){
                ## p1 is first col with sum zero
                p1 <- (1:rowclust)[su==0][1]
                ## p2 is first col with largest sum
                p2 <- (1:rowclust)[su==max(su)][1]
                ## m is the first row with a 1 in column p2
                m <- (1:n)[ppr.ma[,p2]==1][1]
                ## Interchange 1 and 0 in row m between cols p1 and p2:
                ppr.ma[m,p1] <- 1
                ppr.ma[m,p2] <- 0
                su <- apply(ppr.ma,2,sum)
            } 
            
            ## given ppr.ma and ppc.ma update Ym
            su <- apply(ppr.ma,2,sum)
            Ym <- diag(1/su) %*% t(ppr.ma) %*% ys.mat %*% ppc.ma %*% diag(1/sv)  
            
            ## given ys.mat and ppr.ma update ppc.ma
            ppc.ma <- matrix(0,p,colclust)
            Ymu <- ppr.ma %*% Ym
            for (j in 1:p){
                mindif <- sum((ys.mat[,j]-Ymu[,1])^2)
                posmin <- 1
                for (i in 2:colclust){
                    dif <- sum((ys.mat[,j]-Ymu[,i])^2)
                    if (dif < mindif){
                        mindif <- dif
                        posmin <- i
                    } 
                }
                ppc.ma[j,posmin] <- 1
            }
            
            sv <- apply(ppc.ma,2,sum)
            
            ## if a class is empty split the largest class
            while (min(sv)==0){
                ## p1 is first col with sum zero
                p1 <- (1:colclust)[sv==0][1]
                ## p2 is first col with largest sum
                p2 <- (1:colclust)[sv==max(sv)][1]
                ## m is the first row with a 1 in column p2
                m <- (1:p)[ppc.ma[,p2]==1][1]
                ## Interchange 1 and 0 in row m between cols p1 and p2:
                ppc.ma[m,p1] <- 1
                ppc.ma[m,p2] <- 0
                sv <- apply(ppc.ma,2,sum)
            } 
            
            ## given ppr.ma and ppc.ma update Ym
            sv <- apply(ppc.ma,2,sum)
            Ym <- diag(1/su) %*% t(ppr.ma) %*% ys.mat %*% ppc.ma %*% diag(1/sv)
            
            ## Check: is this a better solution?
            B <- ppr.ma %*% Ym %*% t(ppc.ma)
            f <- sum(diag((t(B) %*% B)))/st   
            fdif <- f-f0
            ifelse(fdif > eps,f0 <- f,break)
        }
        if (loop==1){
            Vdkm <- ppc.ma
            Udkm <- ppr.ma
            Ymdkm <- Ym
            fdkm <- f
            loopdkm <- 1
            indkm <- it
            fdifo <- fdif
        }
        if (f > fdkm){
            Vdkm <- ppc.ma
            Udkm <- ppr.ma
            Ymdkm <- Ym
            fdkm <- f
            loopdkm <- loop
            indkm <- it
            fdifo <- fdif
        }
    }
    ## sort clusters of variables per descending order of cardinality
    Vdkm <- Vdkm[,sort.list(-apply(Vdkm,2,sum))]
    ## sort clusters of objects in descending order of cardinality
    Udkm <- Udkm[,sort.list(-apply(Udkm,2,sum))]
    ## Output ppr and ppc:
    list("pp1"=Udkm,"pp2"=Vdkm,"expl.var"=fdkm)
}

