\name{findRQR}
\alias{findRQR}

\title{Find randomised quantile residuals for an already-fitted model}

\description{

  For a given fitted 'clustglm' model with discrete data the uniform
  randomised quantile residuals (URQRs) are chosen using a uniform
  distribution for a distance up the riser in the cumulative
  distribution function (a step function). This essentially jitters the
  value to somewhere in the range of possible values. Re-running the
  'findRQR' function will result in slightly different values because of
  the randomisation. If the model fits the data well, the randomised
  quantile residuals will follow a Uniform[0,1] distribution. The
  equivalent normal randomised quantile residuals (NRQRs) are obtained
  by transforming the URQRs through the inverse standard normal
  distribution function. If the model fits the data well, the NRQRs will
  follow a standard normal distribution.
  
}

\usage{

   findRQR(model)

}

\arguments{

  \item{model}{a model of class 'clustglm' with discrete data
    (e.g. Poisson or binomial family)}
    
}

\details{

  This function is used for a discrete-data model. With continuous data
  use 'findQR', which uses an exact value, the probability integral
  transform. No randomisation is required.
    
}

\value{

  A data frame of unifRQR and normRQR values, one row for each data point.

}

\seealso{

\code{\link{clustglm}}
\code{\link{findQR}}

}

\examples{

## To be provided.
## Reference to Dunn and Smyth to be provided.

}

\keyword{ residuals }

