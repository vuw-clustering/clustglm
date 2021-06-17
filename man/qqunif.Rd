\name{qqunif}
\alias{qqunif}

\title{Uniform Q-Q plot}

\description{

Plot of quantile residuals versus theoretical quantiles from the
Unif[0,1] distribution

}

\usage{

qqunif(y, main = "Uniform Q-Q Plot",
                   xlab = "Theoretical Quantiles",
                   ylab = "Sample Quantiles")

}

\arguments{
  \item{y}{(Randomised) quantile residuals}
  \item{main}{character title (optional)}
  \item{xlab}{character label for x-axis}
  \item{ylab}{character label for y-axis}

}

  
\details{

  'qqunif' plots the (randomised) quantile residuals versus the
  theoretical quantiles from the Uniform[0,1] distribution. For
  continuous data the quantile residual is used, while for discrete data
  a randomised quantile residual is used, with the value being randomly
  chosen from a uniform random variate part-way up the riser in the step
  function. 
  
}

\value{

  The plot is displayed on screen.

}

\seealso{

\code{\link{ordplot}}

}

\examples{

## to be provided
## Reference to Dunn and Smyth to be provided.

}

\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
