\name{logit.fn}
\alias{logit.fn}

\title{Adjusted logit function}

\description{

  This function is like the logit, but replaces arguments too near 0 or
  1, thus bounding the logit away from +/- inf.
  
}

\usage{

logit.fn(x)

}

\arguments{

  \item{x}{any scalar, vector or matrix x}
    
}

\details{

  The value of this function is bounded between +/- 20.92327. This does
  not usually affect the clustering obtained or the log likelihood by
  much. The function is mainly used internally in \code{clustglm}.
  
}

\value{

  A scalar, vector or matrix with values bounded between +/- 20.92327.
  
}

