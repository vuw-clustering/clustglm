## profileplot.R

## Use estimated gamma (interaction) values from the linear predictor.
## Average over any other factors.

profileplot <- function(model, x.factor, trace.factor,
                        main = NULL,
                        sub = NULL,
			xlim = NULL, ylim = NULL,
                        xlab = NULL, ylab = NULL,
                        pch = NULL, lty = NULL, lwd = NULL,
                        xnames = NULL, tracenames = NULL,
                        colourvector = NULL, meanline = FALSE,
			is.main.x = TRUE, sort.x = NULL, ntries = 100000,
                        legend = FALSE, legx = NULL, legy = NULL,
			verbose = 0){
    this.model <- model
    long.df <- this.model$data
    ## Find values of factors:
    x.f <- long.df[,(1:ncol(long.df))[x.factor == names(long.df)]]
    trace.f <- long.df[,(1:ncol(long.df))[trace.factor == names(long.df)]]
    nx <- length(unique(x.f))
    ntr <- length(unique(trace.f))

    ## Find the gamma matrix for plotting:
    these.pars <- findpars(this.model)
    plotpoints.mat <- these.pars[[length(these.pars)]]

    ## R may have reversed factors in the interaction term.
    ## Needs to be nx by ntr:
    if ((nrow(plotpoints.mat) == ntr)&(ncol(plotpoints.mat) == nx))
        plotpoints.mat <- t(plotpoints.mat)

    ## Sort x factor if requested:
    if (!is.null(sort.x)){
        ## Alphanumeric sort of factor names:
        if (sort.x == "alphanum")
	    plotpoints.mat <- plotpoints.mat[sort.list(rownames(plotpoints.mat)),]
        ## Travelling salesperson sort:
        if (sort.x == "TS"){
            ## Initialize the distance travelled (tp = total path) with cols (traces) 1:ntr:
            tp <- 0
            rs <- 1:nx   ## Initial row sort is just 1:nx.
            for(tr in 1:ntr) tp <-
                tp + sum(abs(diff(plotpoints.mat[rs,tr])))
            total.path <- tp    ## Best total path so far
            rows.sort <- rs
            ## Find furthest apart rows:
            dist.mat <- as.matrix(dist(plotpoints.mat,method="manhattan"))
            endpts <- which(dist.mat==max(dist.mat),arr.ind=T)
            start.row <- endpts[1,1]
            end.row   <- endpts[1,2]
            other.rows <- (1:nx)[-c(start.row,end.row)]
            ## Try random orderings of internal columns:
            if (verbose > 0){
		cat("Using travelling salesperson method for ordering on x axis\n")
		cat("Trying ",ntries," random orderings, looking for minimum path\n")
	    }
            for (i in 1:ntries){
	        rs <- c(start.row,sample(other.rows),end.row)
                tp <- 0
                for(tr in 1:ntr) tp <-
                    tp + sum(abs(diff(plotpoints.mat[rs,tr])))
                if (tp < total.path){
                    if (verbose > 0){
   	                cat("New best, total path = ",tp,"\n")
                        cat("Sorted x factor levels:",rs,"\n")
	            }
                    rows.sort <- rs
		    total.path <- tp
                }
            }
            plotpoints.mat <- plotpoints.mat[rows.sort,]
        }
        ## Sorting by one trace:
	if (is.numeric(sort.x)){
	    if (length(sort.x)>1) stop("Can only sort by one level of trace factor")
	    if (length(sort.x)==1)
	        plotpoints.mat <- plotpoints.mat[sort.list(plotpoints.mat[,sort.x]),]
        }
    }   

    ## Set up varying line widths if required:
    if (!is.null(lwd)){
        if (lwd == "vary"){
            if (trace.factor == this.model$clustfactnames[1]){
                pi.v <- this.model$pi.list[[1]]
            }
            if (length(this.model$nclust) == 2){
                if (trace.factor == this.model$clustfactnames[2]){
                    pi.v <- this.model$pi.list[[2]]
                }
            }
            print(pi.v)
            lwd <- pi.v/min(pi.v)  ## Scale 1 to max/min
            print(paste("lwd =",lwd))
        }
    }

    ## Create data frame of plotting points for saving:
    plotpoints.df <- as.data.frame(plotpoints.mat)

    ## Set up for the profile plot:
    this.mat <- plotpoints.mat
    ## Get graphical inputs:
    if (is.null(main)) main = "Profile Plot"
    if (is.null(xlab)) xlab <- x.factor
    if (is.null(ylab)) ylab <- "Interaction coefficient"
    if (is.null(xnames)) xnames <- rownames(this.mat)
    if (is.null(tracenames)) tracenames <- colnames(this.mat)
    if (is.null(colourvector))
        colourvector <- c("red","blue","black","orchid3","brown","green","orange")
    if (is.null(pch)) pch <- rep(16,ntr)
    if (!is.null(pch)&(length(pch) == 1)) pch <- rep(pch,ntr)
    if (!is.null(pch)&(length(pch) >= ntr)) pch <- pch
    ## Warn if 1 < length(pch) < ntr?
    if (is.null(lty)) lty <- rep(1,ntr)
    if (!is.null(lty)&(length(lty) == 1)) lty <- rep(lty,ntr)
    if (!is.null(lty)&(length(lty) >= ntr)) lty <- lty
    ## Warn if 1 < length(lty) < ntr?
    if (is.null(lwd)) this.lwd <- rep(1,ntr)
    if (!is.null(lwd)&(length(lwd) == 1)) this.lwd <- rep(lwd,ntr)
    if (!is.null(lwd)&(length(lwd) >= ntr)) this.lwd <- lwd
    ## Warn if 1 < length(lwd) < ntr?

    ## Plot the profile:
    plot(rep(1:nx,ntr),c(this.mat),type='n',
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, sub = sub,
         axes = F)
    axis(1,1:nx,xnames)
    axis(2)
    box()
    ## For each trace:
    for (tr in 1:ntr){
        points(1:nx, this.mat[,tr], pch = pch[tr], col = colourvector[tr])
        lines(1:nx, this.mat[,tr], lty = lty[tr], lwd = lwd[tr],
              col = colourvector[tr])
    }
    if (meanline == TRUE) lines(c(1,nx), rep(mean(this.mat),2), col = "grey")
    if (legend == TRUE){
        if (is.null(legx)) legx <- round(nx*0.7)
        if (is.null(legy)) legy <- max(this.mat)
        legend(legx,legy,pch=rep(16,ntr),col=colourvector[1:ntr],
               lty=rep(1,ntr),tracenames)
    }
    ## Return the data frame:
    return(plotpoints.df)
}

