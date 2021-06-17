\name{profileplot}
\alias{profileplot}

\title{Profile plot}

\description{

Profile plot of either row groups or column groups of a finite-mixture clustered model.

}

\usage{

profileplot(model, x.factor, trace.factor, main = NULL,
                 sub = NULL, xlim = NULL,
                 ylim = NULL, xlab = NULL, ylab = NULL, pch =
                 NULL, lty = NULL, lwd = NULL, xnames = NULL,
                 tracenames = NULL, colourvector = NULL, meanline =
                 FALSE, is.main.x = TRUE, sort.x = NULL, ntries =
                 100000, legend = FALSE, legx = NULL, legy = NULL,
                 verbose = 0)

}

\arguments{
  \item{model}{clustglm model object model }
  \item{x.factor}{character, name of factor for x axis}
  \item{trace.factor}{character, name of clustered factor for profiles}
  \item{main}{character title (optional)}
  \item{sub}{character subtitle (optional)}
  \item{xlim}{x-axis limits, numeric vector}
  \item{ylim}{y-axis limits, numeric vector}
  \item{xlab}{character label for x-axis}
  \item{ylab}{character label for y-axis}
  \item{pch}{numeric, plotting character}
  \item{lty}{numeric vector, line types for the traces}
  \item{lwd}{numeric vector, line widths for the traces}
  \item{xnames}{character vector, names at ticks along the x axis}
  \item{tracenames}{character vector of cluster names for trace factor}
  \item{colourvector}{character vector of colour names, optional}
  \item{meanline}{logical, if TRUE a horiontal line is drawn at the mean}
  \item{is.main.x}{logical, if FALSE the x factor is not a main effect}
  \item{sort.x}{character "mean" (ascending mean order), "cluster"
  (ascending order for cluster specified), "travellingsalesperson" for
  shortest paths)}
  \item{ntries}{integer, number of random internal orderings to try if
  doing "travellingsalesperson" method for shortest path}
  \item{legend}{logical, if TRUE a legend is plotted}
  \item{legx}{numeric, x value for left edge of legend position}
  \item{legy}{numeric, y value for top of legend position}
  \item{verbose}{numeric; if 0 then intermediate messages and results
  are suppressed, if 1 some output to screen, if 2 copious output to screen}
}

  
\details{

  A profile plot is essentially an interaction plot of the terms in
  A:Bclust or B:Aclust, with any main effects in the model removed. For
  example, in an abundance matrix with 30 sites by 25 species, a
  clustering into say four species clusters yields a 30 by 4 matrix of
  interaction parameters in the linear predictor. This means each of the
  four clusters has a vector of length 30, which may be plotted as an
  index plot, with 1:30 on the x axis and the vector providing the y
  coordinates. On the same plot, there are 30 points on the x axis for
  the levels of the site factor, and four traces or profiles, one for
  each species cluster.
  The sites may be sorted to provide a clearer picture, so that adjacent
  sites have similar species composition. Of the "sort.x" choices, the
  "travellingsalesperson" option puts the least alike sites at the extreme
  ends and tries to find an ordering of the remaining sites which
  overall minimises the total path length of the four paths. The
  travelling salesperson does not have to return to her or his starting
  site. 
  
}

\value{

  The profile plot is displayed on screen and (optionally) the plotting
  points are stored in a specified data frame.

}

\seealso{

\code{\link{ordplot}}

}

\examples{

## to be provided

}

\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
