\name{findstart}
\alias{findstart}

\title{For balanced data, find starting allocations to clusters}

\description{
  
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  For a matrix of data, other clustering methods may be used to quickly
  give initial exact clustering(s) of the rows and/or the columns. This is
  faster than using finite-mixture random starts and often gives useful
  initial allocations before the EM algorithm. Balanced data supplied
  as a data frame is first converted to a matrix.
  
  
}

\usage{

  findstart(x, clustmode = "rows", nclust,
  clusttype = "kmeans", randstarts = 1, metric = "euclidean", dendro = FALSE)

}
%- maybe also 'usage' for other objects documented here.

\arguments{
  
  \item{x}{Either a data matrix for which rows and/or columns are to be
    clustered or a data frame with a response variable and two balanced
    factors for clustering}

  \item{clustmode}{A choice of "rows" or "cols" for single-mode
  clustering, or "both" for biclustering}

  \item{nclust}{numeric vector giving number(s) of clusters (length one
    for single-mode clustering or two for biclustering)}
  
  \item{clusttype}{"kmeans", "hclust" or "diana" for matrix
    row-clustering using k-means (see \code{\link[stats]{kmeans}}),
    hierarchical agglomerative (see \code{\link[stats]{hclust}}) or
    divisive (see \code{\link[cluster]{diana}}) clustering respectively.
    "dkm" does biclustering of rows and columns simultaneously, using
    double k-means}
  
  \item{randstarts}{The number of random starts if \code{kmeans} or
    \code{dkm} is used}
  
  \item{metric}{In \code{hclust} and \code{diana} the default metric
    "euclidean" may be replaced by "manhattan".}
  
  \item{dendro}{If TRUE, the hierarchical methods \code{hclust} and
    \code{diana} will draw a dendrogram.}
  
}

\details{

  These methods have fast clustering algorithms and often provide useful
  initial allocations before deploying the slower, likelihood-based EM
  algorithm. Data should be supplied as a matrix, in which there is one
  observation in each cell(i,j). If a data frame including two balanced
  factors is supplied, conversion to a matrix is done first, with the
  entry in cell(i,j) being the mean of all responses with Factor A at
  level i and B at level j. The number of such responses should be a
  constant, r, the replication number. Once found, these initial
  allocations may be fed into the function \code{clustglm()}.

}

\value{

  A list with one or two allocation-matrix (indicator-matrix)
  components. The rows are the levels of the factor to be clustered and
  the columns are the levels of the associated clustered factor. The
  value in the matrix is 1 if that factor level is in that cluster,
  otherwise zero. With single-mode clustering the list has one component
  giving the allocation for the factor to be clustered. For biclustering
  the list has two items, one allocation matrix for each factor to be
  clustered. 
  
}

\references{
  R Package "cluster", title "Finding Groups in Data", function
  \code{diana}, DIvisive ANAlysis Clustering
  
}

%% \author{
%%  ~~who you are~~
%% }

\note{ For a greater choice of method and metric, the user may take the
  data matrix to any chosen clustering option and obtain a preliminary
  allocation.
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }

\examples{

## Avi Example

data(avi)
x <- as.matrix(avi$fau)
## Row clustering:
pp.start <- findstart(x, clustmode = "rows",
                      nclust = 2,
                      clusttype = "kmeans",
                      randstarts = 100)
## Biclustering:
pp.start <- findstart(x, clustmode = "both",
                      nclust = c(2,2),
                      clusttype = "dkm",
                      randstarts = 100)

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%   \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%   \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
