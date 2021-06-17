\name{randpart}
\alias{randpart}

\title{Random partition of \code{n} objects into \code{c} classes}

\description{

  The number of classes must be less than the number of objects. The
  classes returned will all be non-empty.
  
}

\usage{

   randpart(n,c)

}

\arguments{

  \item{n}{the number of objects}
  \item{c}{the number of classes}
    
}

\details{

  This function is used internally in \code{\link{dkm}} to provide
  starting partitions for double k-means.
  
}

\value{

  An \code{n} by \code{c} indicator matrix, where element [i,j] = 1
  indicates object i is in class j, otherwise zero for
  non-membership. The row sums are all one.

}

\references{

Vichi, M. (2000) Double k-means clustering for simultaneous
classification of objects and variables. In Borra, S. et al. (eds),
\emph{Advances  in  Classification  and  Data  Analysis}, 43-52,
Springer, Berlin. 

}

\seealso{

\code{\link{dkm}}

}


\examples{

library(clustglm)
randpart(12,3)

} 


\keyword{ dkm }

