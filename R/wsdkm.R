## wsdkm.R

## ----------------------------------
## Weighted single or double k-means:
## ----------------------------------

wsdkm <- function(in.df, nAclust, nBclust, randstart, verbose){

    ## Set up notation, objects:
    ## -------------------------
    
    maxiter <- 1
    Qres <- in.df$Qres    ## Normalised (expected) quantile residual
    Nobs <- length(Qres)
    facA <- in.df$facA    ## Given factor (e.g. rows from a matrix)
    facB <- in.df$facB    ## Given factor (e.g. cols from a matrix)
    A <- length(levels(facA))
    B <- length(levels(facB))
    U <- nAclust          ## Number of clusters of facA levels
    V <- nBclust          ## Number of clusters of facB levels
    this.df <- data.frame(Qres, facA, facB)

    ## Do the specified number of random starts:
    ## -----------------------------------------
    
    for (loop in 1:randstart){
    
        ## Generate random partition(s) for the clustered factor(s):
	## ---------------------------------------------------------

	if ((verbose > 0)&(loop == 1))
	    cat("Generating ",randstart,"random starts\n")
	    
	## Factor U (e.g. row clusters), :
	if (U == A) pp1 <- diag(rep(1,A))  ## Trivial case, no clustering
	if (U < A) pp1 <- randpart(A,U)    ## Some clustering
	facU <- rep("u1",Nobs)
	levs.facU <- paste("u",1:U,sep="")
	for (a in 1:A) for (u in 1:U) if (pp1[a,u] == 1)
	    facU[as.numeric(facA) == a] <- levs.facU[u]
	facU <- as.factor(facU)
	this.df$facU <- facU
	
	## Factor V (e.g. column clusters):
	if (V == B) pp2 <- diag(rep(1,B))  ## Trivial case, no clustering
	if (V < B) pp2 <- randpart(B,V)    ## Some clustering
	facV <- rep("v1",Nobs)
	levs.facV <- paste("v",1:V,sep="")
	for (b in 1:B) for (v in 1:V) if (pp2[b,v] == 1)
	    facV[as.numeric(facB) == b] <- levs.facV[v]
	facV <- as.factor(facV)	
	this.df$facV <- facV
	
	## Find centroids and explained variance:
	## --------------------------------------

    	## Construct U by V matrices of cell means, sizes and sums:
    	## Use tapply
	UVmeans <- tapply(X = Qres, INDEX = list(facU,facV),FUN = "mean")
	UVcounts <- tapply(X = Qres, INDEX = list(facU,facV),FUN = "length")
	UVsums <- tapply(X = Qres, INDEX = list(facU,facV),FUN = "sum")
	## NB: UVsums = UVmeans * UVcounts
	## Find fitted for each observation (mean of Qres in cell(u,v)):
	Qres.fits <- rep(NA,Nobs)
	for (uu in 1:U) for (vv in 1:V)
	    Qres.fits[(facU == levels(facU)[uu])&
	              (facV == levels(facV)[vv])] <- UVmeans[uu,vv]
        table(Qres.fits)
	this.df$fits <- round(Qres.fits,3)
	## Total SS:
	this.df$TSSn <- (this.df$Qres - mean(this.df$Qres))^2   ## n-th component
	TSS <- sum(this.df$TSSn)
	## SSE:
	this.df$SSEn <- (this.df$Qres - this.df$fits)^2   ## n-th component
	SSE <- sum(this.df$SSEn)
	## Explained variance:
	EV <- (TSS - SSE)/TSS
	## Starts really low. Hope to improve with iterations.
	EV1 <- round(EV,3)
	bestEV <- EV

        ## Iteration phase:
	## ----------------
	
	## Update clustering by moving points to nearest centroid then update
	## centroids using this clustering.
	## For biclustering, alternate between U and V clusters for compromise.

        it <- 0
        prev.EV <- 0
        this.EV <- EV
        
        while ((it<=maxiter)&(this.EV > prev.EV)){
            it <- it+1
            prev.EV <- prev.EV

	    ## Update A clustering (pp1):
	    ## --------------------------
	    ## Move each level of facA to another cluster if nearer to the centroid.
	    ## Use previous UVmeans matrix and hence previous fitted values.
	    if (U < A){
	        new.pp1 <- pp1
	        ## For each level of A, see if assigning it to a different level of U
		## gives a lower SSE:
		for (a in 1:A){
		    this.SSE <- rep(0,U)
		    small.df <- this.df[this.df$facA == levels(this.df$facA)[a],]
		    for (new.u in 1:U){
		        new.fits <- rep(NA,nrow(small.df))
		        for (v in 1:V){
		            new.fits[small.df$facV == levels(facV)[v]] <-
			        UVmeans[new.u,v]
			}
		        this.SSE[new.u] <- sum((small.df$Qres - new.fits)^2)
	            }
		    new.pp1[a,] <- as.numeric(this.SSE == min(this.SSE))
		}
                ## If a class is empty split the largest class
	    	sum1 <- apply(new.pp1,2,sum)
            	while (min(sum1)==0){
                    ## p1 is first col with sum zero
                    p1 <- (1:U)[sum1==0][1]
                    ## p2 is first col with largest sum
               	     p2 <- (1:U)[sum1==max(sum1)][1]
                     ## m is the first row with a 1 in column p2
                     m <- (1:A)[new.pp1[,p2]==1][1]
                     ## Interchange 1 and 0 in row m between cols p1 and p2:
                     new.pp1[m,p1] <- 1
                     new.pp1[m,p2] <- 0
                     sum1 <- apply(new.pp1,2,sum)
                }
		
		## Given new pp1 and old pp2, update facA and centroids
	    	pp1 <- new.pp1
		facU <- rep("u1",Nobs)
		levs.facU <- paste("u",1:U,sep="")
		for (a in 1:A) for (u in 1:U) if (pp1[a,u] == 1)
		    facU[as.numeric(facA) == a] <- levs.facU[u]
		    facU <- as.factor(facU)
		UVmeans <- tapply(X = Qres, INDEX = list(facU,facV),FUN = "mean")
	    }

	    ## Update B clustering (pp2):
	    ## --------------------------
	    ## Move each level of facB to another cluster if nearer to the centroid.
	    ## Use previous UVmeans matrix and hence previous fitted values.
	    if (V < B){
	        new.pp2 <- pp2
	        ## For each level of B, see if assigning it to a different level of V
		## gives a lower SSE:
		for (b in 1:B){
		    this.SSE <- rep(0,V)
		    small.df <- this.df[this.df$facB == levels(this.df$facB)[b],]
		    for (new.v in 1:V){
		        new.fits <- rep(NA,nrow(small.df))
		        for (u in 1:U){
		            new.fits[small.df$facU == levels(facU)[u]] <-
			        UVmeans[u,new.v]
			}
		        this.SSE[new.v] <- sum((small.df$Qres - new.fits)^2)
	            }
		    new.pp2[b,] <- as.numeric(this.SSE == min(this.SSE))
		}
                ## If a class is empty split the largest class
	    	sum2 <- apply(new.pp2,2,sum)
            	while (min(sum2)==0){
                    ## p1 is first col with sum zero
                    p1 <- (1:V)[sum2==0][1]
                    ## p2 is first col with largest sum
               	     p2 <- (1:V)[sum2==max(sum2)][1]
                     ## m is the first row with a 1 in column p2
                     m <- (1:B)[new.pp2[,p2]==1][1]
                     ## Interchange 1 and 0 in row m between cols p1 and p2:
                     new.pp2[m,p1] <- 1
                     new.pp2[m,p2] <- 0
                     sum2 <- apply(new.pp2,2,sum)
                }
		
		## Given old pp1 and new pp2, update facV and centroids
	    	pp2 <- new.pp2
		facV <- rep("v1",Nobs)
		levs.facV <- paste("v",1:V,sep="")
		for (b in 1:B) for (v in 1:V) if (pp2[b,v] == 1)
		    facV[as.numeric(facB) == b] <- levs.facV[v]
		    facV <- as.factor(facV)
		UVmeans <- tapply(X = Qres, INDEX = list(facU,facV),FUN = "mean")
	    }



	    ## Find SSE and explained variance:
	    ## --------------------------------

	    ## Find fitted for each observation (mean of Qres in cell(u,v)):
	    Qres.fits <- rep(NA,Nobs)
	    for (uu in 1:U) for (vv in 1:V)
	        Qres.fits[(facU == levels(facU)[uu])&
		          (facV == levels(facV)[vv])] <- UVmeans[uu,vv]
            ## table(Qres.fits)
	    this.df$fits <- round(Qres.fits,3)
	    ## Total SS:
	    this.df$TSSn <- (this.df$Qres - mean(this.df$Qres))^2   ## n-th component
	    TSS <- sum(this.df$TSSn)
	    ## SSE:
	    this.df$SSEn <- (this.df$Qres - this.df$fits)^2   ## n-th component
	    SSE <- sum(this.df$SSEn)
	    ## Explained variance:
	    EV <- (TSS - SSE)/TSS
            this.EV <- round(EV,3)
            ## Save if overall best so far:
            if (loop == 1){
                bestEV <- EV
                best.pp1 <- pp1
                best.pp2 <- pp2
            }
            if ((loop > 1)&(EV > bestEV)){
                bestEV <- EV
                best.pp1 <- pp1
                best.pp2 <- pp2
            }
            
            ## Update for next iteration:
            prev.EV <- this.EV
	}
	## Iterations done

	if (verbose > 0){
	    cat("Done randstart number",loop,"\n")
	}
        
    }
    ## End of loops
    if (verbose > 0){
        cat("Done",randstart,"starts\n")
        cat("Best explained variance:",bestEV,"\n")
    
    ## Output non-trivial pp1, pp2,
    if ((U < A)&(V == B)) outlist <- list(pp1 = best.pp1)
    if ((U == A)&(V < B)) outlist <- list(pp2 = best.pp2)
    if ((U < A)&(V < B)) outlist <- list(pp1 = best.pp1, pp2 = best.pp2)
    return(outlist)
}

}
