\name{cottontails}
\docType{data}

\alias{cottontails}
\encoding{UTF-8}
\title{ Cottontails Data }

\description{

  Edwards and Eberhardt (1967) obtained a closed-population
  capture-recapture data set with binary data (1 if captured, 0 if
  uncaptured) for 76 cottontail rabbits (*Sylvilagus floridanus*) over 18
  capture occasions.

  From this, conditioning on the first capture of each rabbit, three
  values were found, "nsucc" = number of recaptures, "ntrials" =
  number of opportunities for recapture and "nfail" = number of failures
  to recapture ("nfail = ntrials - nsucc"). Two rabbits were excluded as
  they were first caught at sample 18 and hence had no opportunity for
  recapture. 

}

\details{

  The recapture data set "cottontails" is a data frame with four
  columns, "ID" = animal identification (a factor), "nsucc" = number of 
  recaptures, "nfail" = number of failures to recapture and "ntrials" =
  number of opportunities for recapture. There are 74 rows, one for each
  animal.  

}

%\source{
%
%}

\seealso{

  \code{\link{clustglm}}

}

\references{

  Edwards, W. R. and Eberhardt, L. (1967) Estimating cottontail
  abundance from livetrapping data.
  *The Journal of Wildlife Management* **31**, 87 -- 96.

}

\examples{


data(cottontails)

ID.2clust.out <-
    clustglm(formula = cbind(nsucc,nfail) ~ IDclust,
             family = "binomial",
             data = cottontails,
             fact4clust = "ID", nclus = 2,
             clustfactnames = "IDclust",
             start.control = list(randstarts = 10),
             verbose = 1)
}

\keyword{datasets}
