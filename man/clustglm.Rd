\name{clustglm}
\alias{clustglm}

\title{Clustering Using Finite Mixtures}

\description{

Likelihood-based GLM models with finite-mixture clustering of factor levels.

}

\usage{

clustglm(formula, family = c("gaussian", "binomial", "poisson"), 
    data, weights = NULL, offset, control = list(...),
    fact4clust = NULL, nclust = NULL, clustfactnames = NULL, 
    start.control = list(randstarts = 0, selfstarts = 0, alloc = NULL),
    EM.control = list(maxEMcycles = 100, EMstoppingpar = 1e-04), 
    save.long = TRUE, save.ests = FALSE, save.Qres = FALSE, save.RQres =
    FALSE, save.EQres = FALSE, verbose = 0, 
    ...)
   } 

\arguments{
  \item{formula}{formula specifying the clustglm model}
  \item{family}{character string naming a family function as in
    \code{\link[stats]{glm}}}
  \item{data}{dataframe (see Details)}
  \item{weights}{prior weights, used for numbers of trials with binomial data}
  \item{offset}{See \code{\link[stats]{glm}}}
  \item{control}{See \code{\link[stats]{glm}}}
  \item{fact4clust}{character vector of existing names for factors to be
    clustered}
  \item{nclust}{numeric vector giving number of clusters (length same as
    \code{clustfactnames})}
  \item{clustfactnames}{character vector of new names for clustered factors}
  \item{start.control}{list of options for initialising the algorithm}
  \item{EM.control}{list of options for controlling the EM algorithm}
  \item{save.long}{if TRUE, save the long data used in the clustglm
    calculations}
  \item{save.ests}{if TRUE, save a sequence of successive
    estimates from the EM algorithm is saved as a matrix \code{ests.seq}}
  \item{save.Qres}{if TRUE, for continuous data e.g. the Gaussian family
    save the quantile residuals (probability integral transforms)}
  \item{save.RQres}{if TRUE, for discrete data e.g. the binomial or
    Poisson families save the randomised quantile residuals}
  \item{save.EQres}{if TRUE, for discrete data (e.g. binomial or Poisson
    families) save the expected quantile residuals, halfway up the riser}
  \item{verbose}{numeric; if 0 then intermediate messages and results
    are suppressed, if 1 some output to screen, if 2 copious output to
    screen}
  \item{\dots}{other arguments}

}

\details{
  The first few arguments, up to \code{control}, are the same as for
  \code{glm()}. The next three arguments specify the clustering:
  \code{fact4clust} names the factor(s) to be clustered, \code{nclust}
  says how many clusters and \code{clustfactnames} defines the new
  clustered factor names. The argument \code{start.control} requires a
  positive number of random starts and/or random self-starts and/or a
  starting allocation to clusters. The argument \code{EM.control} has
  default settings for the maximum number of iterations and the stopping
  rule for convergence of parameter estimates, and so may be incomplete
  or missing.  
 
}

\value{

  \code{clustglm} returns an object of class ``clustglm`` which inherits
  from class ``glm``.  

  The function \code{summary} (i.e., summary.clustglm) can be used to
  obtain or print a summary of the results.

  An object of class ``glm`` is a list containing at least the
  components of a \code{\link[stats]{glm}} object.

  An object of class ``clustglm`` is a list with the components of a
  \code{\link[stats]{glm}} object and a few additional components. These
  comprise both mirrored inputs --

  \item{formula}{formula}
  \item{fact4clust}{fact4clust}
  \item{clustfactnames}{clustfactnames}
  \item{nclust}{nclust}
  \item{start.control}{start.control}
  \item{EM.control}{EM.control}
  
  -- and some other results:
  
  \item{glmformula}{formula stripped of terms that involve clustering}
  \item{pars.record}{if track.pars = TRUE, a matrix tracking the
  successive LL and parameter estimates} 
  \item{pp.list}{a list with matrix components, the posterior
  probabilities of membership for levels of the first (Factor A)
  clustering, and (in the case of biclustering) the posterior
  probabilities of membership for the levels of the second (Factor B)
  clustering. The list is NULL if there is no clustering. With data coming
  from a matrix, Factors A and B may be ROW and COLUMN respectively.}
  \item{final.glm}{the final internal glm model fit in the M step}
  \item{Qres}{a data frame of continuous-data quantile residuals
  (probability integral transforms, PIT)}
  \item{RQres}{a data frame of discrete-data randomised quantile
    residuals}
  \item{EQres}{a data frame of discrete-data expected quantile residuals,
  half-way up the riser}

}

\references{

Pledger, S. and Arnold, R. (2014) Multivariate methods using mixtures:
Correspondence analysis, scaling and pattern-detection. \emph{Computational 
Statistics and Data Analysis} \bold{71}, 241--261.
}

\seealso{
  
  \code{\link{mat2df}}
  \code{\link{findstart}}
  
}

\examples{

## Avi species counts

library(clustglm)
data(avi)
?avi

## Create a long data frame including all the information.

avi.df <- mat2df(y = avi$fau, xr.df = avi$mil, xc.df = avi$traits,
   responsename = "count", factorname1 = "site", factorname2 = "species")
str(avi.df)
head(avi.df)

## Fit a pattern-detection model with two species clusters:

species2clust.out <-
    clustglm(formula = count ~ site + species + site:speciesclust,
             family = "poisson",
             data = avi.df,
             fact4clust = "species", nclust = 2,
             clustfactnames = "speciesclust",
             start.control = list(selfstarts = 10),
             verbose = 1)

## Allocate each species to a cluster using posterior probabilities:

pp.mat <- species2clust.out$pp.list[[1]]
split(rownames(pp.mat),apply(pp.mat,1,which.max))

}

\keyword{ models }
\keyword{ cluster }

