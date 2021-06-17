## E_step.R

## Used to update posterior probabilities.
## ---------------------------------------

## Can do unbalanced data.

E_step <- function(nf4c, nlev, nclus,
                   f4c.colno, cf.colno,
                   pp.list, pi.list,
                   long.df, verbose){

    if (verbose > 1){
        cat("Starting E step\n")
        cat("Input pp.list\n")
        print(pp.list)
    }

    ## Define constants:
    A <- nlev[1]
    U <- nclus[1]
    if (nf4c == 2){
        B <- nlev[2]
        V <- nclus[2]
    }
    
    ## For within E step, extend long data frame with factors A,U,B,V:

    E.df <- data.frame(long.df,
                       A = long.df[,f4c.colno[1]],
                       U = long.df[,cf.colno[1]])
    if (nf4c == 2){
        ## For biclustering, 
        E.df$B <- long.df[,f4c.colno[2]]
        E.df$V <- long.df[,cf.colno[2]]
    }

    ## Single-mode clustering:
    ## -----------------------
    
    if (nf4c == 1){
        cf <- 1
        ## If non-trivial:
        if ((nclus[cf]>1)&(nclus[cf]<nlev[cf])){
            for (a in 1:nlev[cf]){
                A.df <- E.df[which(E.df$A == levels(E.df$A)[a]),]
                lognum <- rep(NA,nclus[cf])
                for (u in 1:nclus[cf]){
                    AU.df <- A.df[which(A.df$U == levels(A.df$U)[u]),]
                    lognum[u] <- log.fn(pi.list[[cf]][u]) +
                        sum(AU.df$logpdf)
                }
                adj.lognum <- lognum - max(lognum,na.rm=T)
                num <- exp(adj.lognum)
                denom <- sum(num,na.rm=T)
                pp.list[[cf]][a,] <- num/denom
            }
            pp.list[[cf]][is.na(pp.list[[cf]])] <- 0
        }
    }        
 
    ## Biclustering:
    ## -------------
    
    ## For each clustered factor in turn, update pp.list.
    
    if (nf4c == 2){
        ## First clustered factor:
        if (verbose > 1)
            cat("Biclustering: updating first clustered factor\n")
        cf <- 1
        ## If non-trivial:
        if ((nclus[cf]>1)&(nclus[cf]<nlev[cf])){
            pp1.mat.new <- matrix(NA,nlev[cf],nclus[cf])
            for (a in 1:nlev[cf]){
                A.df <- E.df[which(E.df$A == levels(E.df$A)[a]),]
                A.df$y.index <- droplevels(A.df$y.index)
                lognum <- log.fn(pi.list[[1]])
                for (u in 1:U){
                    AU.df <- A.df[which(A.df$U == levels(A.df$U)[u]),]
		    ## Will have same yi values as A.df.
		    lognum[u] <-lognum[u] +
		        sum(log.fn(tapply(AU.df$pp2.pdf, AU.df$y.index, sum)))
		}
                adj.lognum <- lognum - max(lognum,na.rm=T)
                num <- exp(adj.lognum)
                denom <- sum(num,na.rm=T)
                pp1.mat.new[a,] <- num/denom
            }
            
            pp1.mat.new[is.na(pp1.mat.new)] <- 0
        }

        if (verbose > 1){
            temp <- round(pp1.mat.new,3)
	    cat("Output pp1\n")
            print(temp)
        }

        ## Second clustered factor:
        if (verbose > 1)
            cat("Biclustering: updating second clustered factor\n")
        cf <- 2
        ## If non-trivial:
        if ((nclus[cf]>1)&(nclus[cf]<nlev[cf])){
            pp2.mat.new <- matrix(NA,nlev[cf],nclus[cf])
            for (b in 1:nlev[cf]){
                B.df <- E.df[which(E.df$B == levels(E.df$B)[b]),]
                B.df$y.index <- droplevels(B.df$y.index)
                lognum <- log.fn(pi.list[[2]])
                for (v in 1:V){
                    BV.df <- B.df[which(B.df$V == levels(B.df$V)[v]),]
		    ## Will have same yi values as B.df.
		    lognum[v] <-lognum[v] +
		        sum(log.fn(tapply(BV.df$pp1.pdf, BV.df$y.index, sum)))
		}
                adj.lognum <- lognum - max(lognum,na.rm=T)
                num <- exp(adj.lognum)
                denom <- sum(num,na.rm=T)
                pp2.mat.new[b,] <- num/denom
            }
            
            pp2.mat.new[is.na(pp2.mat.new)] <- 0
        }
	
        if (verbose >1){
            temp <- round(pp2.mat.new,3)
	    print("Output pp2")
            print(temp)
        }
	
        pp.list[[1]] <- pp1.mat.new
        pp.list[[2]] <- pp2.mat.new
    }           
    ## Finished biclustering section.
    
    if (verbose > 1) cat("Finished E step\n")  
    return(pp.list)
}
