\name{ordplot}
\alias{ordplot}

\title{

  Ordination plot

}

\description{

  Create an ordination plot, using the clustering of one factor to
  ordinate levels of the other factor (e.g. cluster rows to ordinate
  columns).
  

}

\usage{

ordplot(model, ord.factor, clust.factor, main = NULL, sub = NULL,
        xlim = NULL, ylim = NULL, xlab = "SVD1", ylab = "SVD2",
        pch = 4, pt.labels = NULL, colourvector = NULL)
}

\arguments{
  \item{model}{clustglm model object model }
  \item{ord.factor}{character}
  \item{clust.factor}{character}
  \item{main}{character title (optional)}
  \item{sub}{character subtitle (optional)}
  \item{xlim}{x-axis limits, numeric vector}
  \item{ylim}{y-axis limits, numeric vector}
  \item{xlab}{character label for x-axis}
  \item{ylab}{character label for y-axis}
  \item{pch}{plotting character for points}
  \item{pt.labels}{optional character vector}
  \item{colourvector}{character vector of colour names, optional}
}

\details{

  The ordination plot uses the interaction terms of the model to display
  closeness of the levels of the unclustered factor (ordination factor,
  ord.factor) using the clustered factor to provide coordinates. Each
  level of the ordination factor is plotted as a point. In general, two
  clusters of the clustered factor provides ordination of the full
  factor in two dimensions. However for the family "poisson" there is a
  loss of one dimension, so three clusters gives ordination in 3D but
  the points are coplanar and so the plane is rotated down to a 2D
  plot. For Poisson data, use ordplot() on a model with three clusters
  and for other families use two clusters.
  

}

\value{

  The ordination plot is displayed on screen and (optionally) the
  plotting points are stored in a specified data frame.

}

\seealso{

\code{\link{clustbiplot}}

}

\examples{


}

