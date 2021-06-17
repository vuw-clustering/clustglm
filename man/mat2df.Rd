\name{mat2df}
\alias{mat2df}

\title{Convert matrix to long data frame}

\description{

  A data matrix or a matrix-style data frame of response values is
  converted to the long dataframe format used for input to
  \code{\link{clustglm}}. The responses are reformatted from a matrix to
  a single column in a long data frame. Optionally, data frames of row-based
  covariates and/or column-based covariates may be included: these will
  provide extra columns in the output data frame. For binomial data the
  response is the number of successes but an extra matrix or
  matrix-style data frame of number of trials is required. Its
  dimensions must match the response input.
  
}

\usage{

  mat2df(y, ntrials = NULL, xr.df = NULL, xc.df = NULL,
  responsename = "Y", factorname1 = "facA", factorname2 = "facB")

}

\arguments{

  \item{y}{numeric data matrix or data frame; data must all be of
    the same type, e.g. all binary, all counts, all continuous}
  \item{ntrials}{if the y data are numbers of successes from binomial
    distributions, a numeric data matrix or data frame matching y and
    giving the number of trials for each count in y}
  \item{xr.df}{optional dataframe of row covariates; row names must
    match row names of y, columns may be nominal, discrete or
    continuous covariates}
  \item{xc.df}{optional dataframe of column covariates; row names must
    match column names of y, columns may be nominal, discrete or
    continuous covariates}
  \item{responsename}{optional name of the intended response variable;
    default is Y}
  \item{factorname1}{optional name of the y row factor in the output data
    frame; default is "facA"}
  \item{factorname2}{optional name of the y column factor in the
    output data frame; default is "facB"}
}

\details{

  The use of "glm" requires the data to be in a data frame with one row
  per observed response. This function stretches an A by B matrix of
  responses into a data frame with AB rows.

}

\value{

  Dataframe with at least 3 columns ("Y", "facA", "facB", or optional
  supplied names), and one additional column for each covariate from 
  \code{xr.df} and \code{xc.df}. For binomial responses, another two
  columns are 'ntrials' and 'nfails'.

}

\seealso{

\code{\link{clustglm}}

}

\examples{
## Count data, avi species counts

data(avi)
?avi
avi.df <- mat2df(y = avi$fau, xr.df = avi$mil, xc.df = avi$traits,
   responsename = "count", factorname1 = "site", factorname2 = "species")
str(avi.df)
head(avi.df)

## Binomial data, avi binary site features

data(avibinary)
?avibinary
avibin.df <- mat2df(y = avibinary,
   responsename = "count", factorname1 = "site", factorname2 = "species")


## If instead the data are from binomial distributions, the a by b
## response matrix is the number of successes but a matrix of the same
## size for the number of trials must also be provided. The exception is
## for binary 0/1 data where (with no ntrials matrix given) the number
## of trials is assumed to be always one.


}

\keyword{ manip }

