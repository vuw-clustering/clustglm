\name{M_step}
\alias{M_step}

\title{M step for clustglm}

\description{

  Mostly hidden from view and used inside \code{clustglm}, this function
  uses the long data frame and other information to run the M step in
  the EM algorithm. It updates the parameter estimates using the latest
  posterior probability estimates from the E_step, and uses \code{glm}
  to fit the model assuming the current complete data estimates. 
}

\usage{

M_step(family, formula, response, nf4c, nlev, nclus, f4c.colno,
cf.colno, indcf.list, pp.list, long.df, verbose)

   } 

\arguments{

  \item{family}{character string naming a family function as in
    \code{\link[stats]{glm}}}
  \item{formula}{\code{clustglm} formula}
  \item{response}{response in the formula, usually a column of the data
    frame but with binomial data may be a matrix of two columns, numbers
    of successes and failures}
  \item{nf4c}{number of factors for clustering}
  \item{nlev}{number(s) of levels of factor(s) to be clustered}
  \item{nclus}{number(s) of levels of clustered factor(s)}
  \item{f4c.colno}{column number(s) in long data frame of the factor(s)
    for clustering}
  \item{cf.colno}{vector of column number(s) in long data frame of the
    clustered factor(s)}
  \item{indcf.list}{list of indicator matrices}
  \item{pp.list}{list of current posterior probability matrices}
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

  \code{M_step} returns a list of (i) pi.list = list of estimated
  probabilities of levels of each factor-fo-clustering being allocated
  to each level of the clustered factor, (ii) pp.list = current estimate
  of posterior probabilities, (iii) LLc =  current estimate of log
  likelihood under complete knowledge, (iv) LLint = current estimate of
  integrated (marginal) log likelihood, (v) long.df = long data frame
  holding all the original data and current allocations, and (vi)
  this.glm = latest \code{glm()} fit assuming the current complete data
  estimates.

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
  \code{\link{E_step}}
  
}
