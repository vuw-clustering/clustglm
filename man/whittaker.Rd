\name{whittaker}
\docType{data}

\alias{whittaker}

\title{ Whittaker Shrub Data }

\description{

Counts of shrubs and seedlings of 25 species over 10 sites in the Siskiyou mountains, Oregon and California.

}

\details{

Whittaker (1960, p.289) provided counts of \emph{J} = 25 shrub and seedling
species over \emph{I} = 10 sites in the Siskiyou Mountains, Oregon and California.

This data set is a list of three objects:

1. whittaker[[1]], "counts", is a matrix of counts of 25 species by 10 sites.
2. whittaker[[2]], "sitecov", is a data frame giving the site covariate
"terrain" with levels "sheltered" and "open".
3. whittaker[[3]], "speciescov", is a data frame giving the species
covariate "Plant.Type" with levels "conifer seedling",
"broadleaf seedling" and "shrub".

The sites 1--10 were chosen by Whittaker in a transect from a valley
floor to a ridge, giving a progression from moist, sheltered sites to
dry, exposed sites. The biggest change observed by Whittaker was between
sites 3 and 4, so we have constructed a site covariate `terrain'
("sheltered" for sites 1--3 and "open" for sites 4--10).

}

\source{

Whittaker (1960, p. 289 Table 4).

}

\seealso{

  \code{\link{clustglm}}

}

\references{

Whittaker, R. H. (1960) Vegetation of the Siskiyou Mountains, Oregon and California. \emph{Ecological Monographs} \bold{30}, 279--338.

}

\examples{

data(whittaker)

## Construct the data frame for \code{clustglm}:

whit.df <- mat2df(y = whittaker$counts,
                  xr.df = whittaker$speciescov,
                  xc.df = whittaker$sitecov,
                  responsename = "count",
                  factorname1 = "species",
                  factorname2 = "site")
head(whit.df)

## Fit a model with species (row) clustering:

whit.clust <- clustglm(formula = count ~ species + site + speciesClust:site,
                       family = "poisson", data = whit.df,
                       fact4clust = "species", nclus = 2,
                       clustfactnames=c("speciesClust"),
                       start.control = list(randstarts = 10),
                       EM.control = list(startEMcycles = 3, EMcycles = 100,
                           EMstoppingpar = 1e-04),
                       verbose = 1)

}

\keyword{datasets}
