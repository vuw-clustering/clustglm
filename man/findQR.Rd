\name{findQR}
\alias{findQR}

\title{Find quantile residuals for an already-fitted model}

\description{

  For a given fitted model with continuous data and the Gaussian family,
  the uniform quantile residuals (UQRs) are the values of F(y), where
  F(.) is the distribution function for y, using the estimated
  parameters. The UQRs are also known as the probability integral
  transforms. If a finite mixture model was fitted, conditional UQRs are
  calculated first using the possible fitted parameters. For each data
  point these conditional UQRs are multiplied by the posterior weights
  for cluster membership, then summed to give a single QR for each
  datum. If the model fits the data well, the quantile residuals will
  follow a Uniform[0,1] distribution. The equivalent normal quantile
  residuals (NQRs) are obtained by transforming the UQRs through the
  inverse standard normal distribution function. If the model fits the
  data well, the NQRs will follow a standard normal distribution.  
  
}

\usage{

   findQR(model)

}

\arguments{

  \item{model}{a model of class 'clustglm' with continuous data
    (e.g. the gaussian family)}
    
}

\details{

  This function is used for a continuous-data model. With discrete data
  use 'findRQR', in which a continuity approximation is done by
  jittering values up the risers of the CDF step functions. 
    
}

\value{

  A data frame of unifQR and normQR values, one row for each data point.

}

\seealso{

\code{\link{clustglm}}
\code{\link{findRQR}}

}

\examples{

## To be provided.
## Reference to Dunn and Smyth to be provided.

}

\keyword{ residuals }

