\name{clustbiplot}
\alias{clustbiplot}

\title{

  Clustered biplot

}

\description{

  Create a biplot with clustering information included. Only implemented
  for family "poisson" so far.

}

\usage{

clustbiplot(model, data.mat, verbose = 1, main = NULL, sub = NULL,
        xlim = NULL, ylim = NULL, xlab = "Dimension 1", ylab =
        "Dimension 2", pch = NULL, rnames = NULL, cnames = NULL,
        rcnames = NULL, ccnames = NULL, colrvect = NULL, myratio = 1,
        centroids = TRUE, partitions = TRUE, axis.lines = TRUE)
}

\arguments{
  \item{model}{clustglm model object with 3 by 3 clustering, family "poisson" }
  \item{data.mat}{matrix of counts}
  \item{verbose}{value 1 gives some output to screen}
  \item{main}{character title (optional)}
  \item{sub}{character subtitle (optional)}
  \item{xlim}{x-axis limits, numeric vector}
  \item{ylim}{y-axis limits, numeric vector}
  \item{xlab}{character label for x-axis}
  \item{ylab}{character label for y-axis}
  \item{pch}{plotting character for points}
  \item{rnames}{optional character vector for row names}
  \item{cnames}{optional character vector for column names}
  \item{rcnames}{optional character vector for row cluster names}
  \item{ccnames}{optional character vector for column cluster names}
  \item{colrvect}{character vector of colour names, optional}
  \item{myratio}{ratio of influence of rows versus columns in the plot}
  \item{centroids}{whether to plot rowcluster and colcluster points}
  \item{partitions}{not yet implemented}
  \item{axis.lines}{grey dashed lines for axis location}
}

\details{

  See Pledger and Arnold 2014. More details to be added here.

}

\value{

  The biplot is displayed on screen and (optionally) the plotting points
  are stored in a specified data frame.

}


\examples{

## to be provided

}

