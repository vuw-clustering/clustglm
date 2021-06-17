\name{counts}
\docType{data}

\alias{counts}

\title{ Count data }

\description{

  Pledger and Arnold (2014) illustrated their mixture-based multivariate
  methods using an 8 by 10 data matrix of counts, with rows labelled A
  to H and columns labelled a to j.

}

\details{

  The counts range from 0 to 16. Patterns of similar rows and columns
  within the matrix are not immediately apparent, but analysis using
  clustglm models shows the block-diagonal structure.
  
}

\seealso{

  \code{\link{clustglm}}

}

\references{

  Pledger, S. and Arnold, R. (2014) Multivariate methods using mixtures:
  Correspondence analysis, scaling and pattern-detection. 
  \emph{Computational Statistics and Data Analysis} \bold{71}, 241--261.

}

\examples{

data(counts)

counts.df <- mat2df(y = counts, responsename = "count",
                    factorname1 = "facA", factorname2 = "facB")
str(counts.df)

## Two facA clusters (row clustering in the original matrix):
facA2clust.out <- clustglm(formula = count ~ facA + facB + facAclust:facB,
                           family = "poisson", data = counts.df,
                           fact4clust = "facA", nclus = 2,
                           clustfactnames = "facAclust",
                           start.control = list(randstarts = 50))
names(facA2clust.out)       
pp.mat <- facA2clust.out$pp.list[[1]]
round(pp.mat,3)
split(rownames(pp.mat), apply(pp.mat,1,which.max))

profileplot(model = facA2clust.out, x.factor = "facB",
            trace.factor = "facAclust",
            main = "Profile plot, two facA clusters",
            sort.x = 1, colourvector = c("red","blue"),
            legend = TRUE, legx = 5, legy = 1)

## Three A clusters:
facA3clust.out <- clustglm(formula = count ~ facA + facB + facAclust:facB,
                         family = "poisson", data = counts.df,
                         fact4clust = "facA", nclus = 3,
                         clustfactnames = "facAclust",
                         start.control = list(randstarts = 50))
names(facA3clust.out)       
pp.mat <- facA3clust.out$pp.list[[1]]
round(pp.mat,3)
split(rownames(pp.mat), apply(pp.mat,1,which.max))

ordplot(model = facA3clust.out, ord.factor = "facB",
        clust.factor = "facAclust", xlab = "Dim1", ylab = "Dim2")

}



\keyword{datasets}
