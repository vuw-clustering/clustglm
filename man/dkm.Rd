\name{dkm}
\alias{dkm}

\title{Double k-means}

\description{

  A matrix is biclustered, i.e. row clusters and column clusters are
  formed simultaneously. The number of random starts is specified and
  for each a least squares method similar to k-means is used. The
  biclustering with the lowest sum of squares is chosen.
  
}

\usage{

   dkm(y.matrix, rowclust, colclust, randstart)

}

\arguments{

  \item{y.matrix}{the numerical matrix to be clustered}
  \item{rowclust}{the number of row clusters required}
  \item{colclust}{the number of column clusters required}
  \item{randstart}{the number of random starts to try}
    
}

\details{

  Each random start begins with a random double partition of the rows
  into row groups and the columns into column groups, using the function
  \code{\link{randpart}}. Next a few iterations of an alternating least squares
  algorithm gives convergence to a better biclustering. The biclustering
  with the minimum sum of squares is chosen. Although this is designed for
  continuous data, it may also provide a useful starting point for
  fitting biclustered models to discrete numerical data.
  
}

\value{

  A list with components \code{pp1} = posterior probabilities for rows into
  row clusters, \code{pp2} for columns into column clusters, and
  code{expl.var} = explained variance.  

}

\references{

Vichi, M. (2000) Double k-means clustering for simultaneous
classification of objects and variables. In Borra, S. et al. (eds),
\emph{Advances  in  Classification  and  Data  Analysis}, 43-52,
Springer, Berlin. 

}

\seealso{

\code{\link{randpart}}

}


\examples{

library(clustglm)
my.mat <- matrix(1:28,4,7)
my.dkm <- dkm(my.mat,rowclust = 2, colclust = 3, randstart = 100)
my.dkm

} 


\keyword{ dkm }

