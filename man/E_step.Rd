\name{E_step}
\alias{E_step}

\title{E step for clustglm}

\description{

  Mostly hidden from view and used inside \code{clustglm}, this function
  uses the long data frame and other information to run the E
  (expectation) step in the EM algorithm. It updates the estimates of
  posterior probabilities using parameters from the \code{glm()} fitted in
  the latest M step.
  
}

\usage{

E_step(nf4c, nlev, nclus, f4c.colno, cf.colno, pp.list, pi.list,
long.df, verbose)

}

\arguments{

  \item{nf4c}{number of factors for clustering}
  \item{nlev}{number(s) of levels of factor(s) to be clustered}
  \item{nclus}{number(s) of levels of clustered factor(s)}
  \item{f4c.colno}{column number(s) in long data frame of the factor(s)
    for clustering}
  \item{cf.colno}{vector of column number(s) in long data frame of the
    clustered factor(s)}
  \item{pp.list}{list of current posterior probability matrices}
  \item{pi.list}{list of current allocation probability estimates}
  \item{long.df}{long data frame with multiple copies of the original
    data, one copy for each combination of factor levels of the
    clustered factors, with extra columns for the current values of
    weights, posterior probabilities}
  \item{verbose}{numeric; if 0 then intermediate messages and results
    are suppressed, if 1 some output to screen, if 2 copious output to
    screen}

}

\details{
  This function is usually used only internally in \code{clustglm},
  where the arguments are supplied internally. However, the information
  above should help if the user looks at the internal details. 

}

\value{

  \code{E_step} returns pp.list, a list of current estimates of
  posterior probabilities.  

}


\references{

  Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977). Maximum
  Likelihood for Imcomplete Data via the EM Algorithm. \emph{Journal of
  the Royal Statistical Society Series B (Methodological)} \bold{39(1)}, 1--38.

  Pledger, S. and Arnold, R. (2014) Multivariate methods using mixtures:
  Correspondence analysis, scaling and pattern-detection. \emph{Computational 
    Statistics and Data Analysis} \bold{71}, 241--261. 
}

\seealso{
  
  \code{\link{clustglm}}
  \code{\link{M_step}}
  
}
