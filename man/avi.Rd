\name{avi}
\docType{data}

\alias{avi}
\encoding{UTF-8}
\title{ Avi Data }

\description{

  A random subset of the "aviurba" data set in the package
  "ade4"

  The full data set "aviurba" has a list including three data frames:

  * "aviurba$fau" is a data frame of abundance (counts) of 40 bird species
  over 51 sites (relevés) near Lyons, France,

  * "aviurba$mil" has eleven environmental variables for the 51 sites,

  * "aviurba$traits" has four traits for the 40 species.

  For the data set "avi" a random sample of 12 species and 15 sites was
  taken, and the "aviurba" data frames reduced accordingly. One
  environmental variable "veg.cover" has been discarded as the
  random sample of sites reduced the number of observed levels of this
  variable. Also, the species were renamed using their abbreviated
  Latin names, e.g. "Aa" for *Alauda arvensis* (skylark) instead of
  "Sp9".
 
}

\details{

  The data set "avi" is a list with four items:

  * "avi$fau" is a data frame of abundance (counts) of 12 bird species
  over 15 sites (relevés, **R1**, **R3**, etc.),

  * "avi$mil" has ten environmental variables for the 15 sites,

  * "avi$traits" has four traits for the 12 species,

  * "avi$species.names" is a matrix with the full names and the name
  codes for the 12 species.
  
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


data(avi)
avi.df <- mat2df(y = avi$fau, xr.df = avi$mil, xc.df = avi$traits,
                     responsename = "count",
                     factorname1 = "site", factorname2 = "species")


species.clust.out <-
    clustglm(formula = count ~ site + species + site:speciesClust,
             family = "poisson",
             data = avi.df,
             fact4clust = "species", nclus = 2,
             clustfactnames = "speciesClust",
             start.control = list(randstarts = 10),
             verbose = 1)

}

\keyword{datasets}
