\name{liphook}
\docType{data}

\alias{liphook}

\title{ Liphook Forest Toadstool Data }

\description{

In an experiment in the Liphook Forest, U.K., Peter Shaw counted ten
species of fungi in seven fixed quadrats over five years, 1986--1990
(Shaw, 2003; Shaw et al., 2003). The data were pooled over (similar)
quadrats by Pledger and Arnold (2014).

}

\details{

The data set is a 10 by 5 matrix of pooled counts of ten species (Bf,
Cs, Gr, Il, Lp, Lr, Pi, Sb, Sl and Sv, initials from their Latin names)
over five years (86, 87, 88, 89, 90 for 1986-1990). The counts range
from 0 to 5823, with only two species present in 1986.
  
}

\source{

  Shaw \emph{et al.} 2003.

}

\seealso{

  \code{\link{clustglm}}

}

\references{

Pledger, S. and Arnold, R. (2014) Multivariate methods using mixtures:
Correspondence analysis, scaling and pattern-detection. 
\emph{Computational Statistics and Data Analysis} \bold{71}, 241--261.

Shaw, P.J.A. (2003) \emph{Multivariate Statistics for the Environmental Sciences}. Arnold, London.

Shaw, P.J.A., Kibby, G. and Mayes, J. (2003) Effects of thinning treatment on an ectomycorrhizal succession under Scots pine. \emph{Mycological Research} \bold{107}, 317--328.

}

\examples{

data(liphook)

liphook.df <- mat2df(y = liphook, responsename = "count",
                     factorname1 = "species",
                     factorname2 = "year")

str(liphook.df)

sp2clust.out <- clustglm(formula = count ~ species + year + speciesClust:year,
                         family = "poisson", data = liphook.df,
                         fact4clust = "species", nclus = 2,
                         clustfactnames = "speciesClust",
                         start.control = list(randstarts = 50),
                         EM.control = list(startEMcycles = 3,
                             EMcycles = 100, EMstoppingpar = 1e-04),
                         verbose = 0)
names(sp2clust.out)       
round(sp2clust.out$pp.list[[1]],3)

}

\keyword{datasets}
