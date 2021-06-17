\name{comparison}
\alias{comparison}

\title{Compare models using AIC, BIC}

\description{

  A list of fitted models is used to form a data frame showing log
  likelihood, number of parameters, AIC, relative AIC, BIC and relative
  BIC. It may be used for model selection and for plotting models with x
  = npar and y = log likelihood. 
  
}

\usage{

   comparison(modellist,rounding = 3)

}

\arguments{

  \item{modellist}{a list of clustglm models applied to the same data set}
    \item{rounding}{The level of rounding for the table}
    
}

\details{

  Model selection by the lowest AIC or lowest BIC is easily obtained
  from the table.
  
}

\value{

  A data frame with the appearance of a table of models and their
  npar, estimated log likelihood, AIC and BIC. This is shown on screen,
  and is stored if assigned to an object.
  

}

\seealso{

\code{\link{clustglm}}

}

\examples{

## To be provided.

}

\keyword{ compare }

