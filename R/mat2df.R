## mat2df.R

mat2df <- function(y, ntrials = NULL, xr.df = NULL, xc.df = NULL,
                   responsename = "Y",
                   factorname1 = "facA",
                   factorname2 = "facB"){
    ## Used to stretch out an A*B data frame or matrix y and up to two
    ## covariate data frames into a data frame with AB observations.
    ## If binomial data, ntrials is a matrix or data frame matching
    ## y with the number of trials in each case.
    ## xr.df is a data frame of row-based covariates (eg. site factors)
    ## xc.df is a data frame of column-based covariates (eg. species traits)
    ## Lengthen the row covariates data frame:
  
  y <- as.matrix(y)  ## coerce to matrix, in case of dataframe

  if (!(is.data.frame(xr.df))&(!is.null(xr.df)))
      stop("Please supply the row covariates as a data frame")
      
  if(is.data.frame(xr.df)){
    xr.df2 <- xr.df
    for (i in 1:(ncol(y)-1))
      xr.df2 <- rbind(xr.df2,xr.df)
  }

  if (!(is.data.frame(xc.df))&(!is.null(xc.df)))
      stop("Please supply the column covariates as a data frame")
      
  ## For xc, construct in wrong order, then rearrange:
  if(is.data.frame(xc.df)){
    xc.df2 <- xc.df
    for (i in 1:(nrow(y)-1))
      xc.df2 <- rbind(xc.df2,xc.df)
    ## Rearrangement preserves factors, numeric
    for (j in 1:ncol(xc.df2))
      xc.df2[,j] <- rep(xc.df[,j],each=nrow(y))
  }
  ## Build new data frame:
  if (is.null(rownames(y))) rownames(y) <- as.character(1:nrow(y))
  if (is.null(colnames(y))) colnames(y) <- as.character(1:ncol(y))
  my.df <- data.frame(Y = c(y),
                      ROW = gl(nrow(y),1,nrow(y)*ncol(y),
                          labels = rownames(y)),
                      COL = gl(ncol(y),nrow(y),nrow(y)*ncol(y),
                          labels = colnames(y)))
  names(my.df)[1] <- responsename
  names(my.df)[2] <- factorname1
  names(my.df)[3] <- factorname2
  if(is.data.frame(xr.df))
    my.df <- cbind(my.df,xr.df2)
  if(is.data.frame(xc.df))
    my.df <- cbind(my.df,xc.df2)
  rownames(my.df) <- as.character(1:nrow(my.df))
  ## If binomial data:
  if (!is.null(ntrials)){
      if (is.data.frame(ntrials)) ntrials <- as.matrix(ntrials)
      if ((nrow(ntrials) != nrow(y))|(ncol(ntrials) != ncol(y)))
          stop("ntrials must be entered as a data frame or matrix matching the dimensions of y")
      my.df$ntrials <- c(ntrials)
      my.df$nfail <- my.df$ntrials - my.df[,1]
  }
  ## If all responses 0 or 1, temporarily assume binomial data:
  if (all(my.df[,1] %in% c(0,1))){
      my.df$ntrials <- rep(1,nrow(my.df))
      my.df$nfail <- my.df$ntrials - my.df[,1]
  }
  ## Return data frame:
  return(my.df)
}
