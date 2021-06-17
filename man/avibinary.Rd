\name{avibinary}
\docType{data}

\alias{avibinary}
\encoding{UTF-8}
\title{Aviurba Binary Data}

\description{

  Sites and site factors from the "aviurba" data set in the package
  "ade4".

  The full data set "aviurba" has a list including three data frames;

  * "aviurba$fau" is a data frame of abundance (counts) of 40 bird species
  over 51 sites (relevés) near Lyons, France,

  * "aviurba$mil" has eleven environmental variables for the 51 sites,

  * "aviurba$traits" has four traits for the 40 species.

  The data set "avibinary" is the data frame "aviurba$mil" reduced to the
  ten binary environmental variables, such as the presence/absence of
  farms, high buildings, etc. The (continuous) percentage vegetation cover
  variable has been removed. 
 
}

\details{

  The data set "avibinary" is a data frame of 51 "sites" (relevés) R1
  ... R51 each with ten binary "features" as detailed in "ade4". These
  habitat attributes describe the degree of urbanization: the presence
  or absence of farms, small buildings, high buildings, industry, fields,
  grassland, scrubby areas, deciduous woods, coniferous woods and noisy
  areas. 
 
  The sites may be clustered using their features, and vice-versa.
  
}

\source{

R package **ade4**

}

\seealso{

  \code{\link{clustglm}}

}

\references{

  Dolédec, S., Chessel, D., Ter Braak, C. J. F. and Champely S. (1996).
  Matching species traits to environmental variables: a new three-table
  ordination method. 
  *Environmental and Ecological Statistics* **3**, 143--166.

}

\examples{

data(avibinary)

avibin.df <- mat2df(y = avibinary,
                 responsename = "pres.abs",
                 factorname1 = "site", factorname2 = "feature")

site2clust.out <-
    clustglm(formula = cbind(pres.abs, ntrials) ~ site + feature + siteclust:feature,
             family = "binomial",
             data = avibin.df,
             fact4clust = "site", nclus = 2,
             clustfactnames = "siteclust",
             start.control = list(randstarts = 10),
             verbose = 1)
}

\keyword{datasets}
