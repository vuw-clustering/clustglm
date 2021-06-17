\name{xlogx.fn}
\alias{xlogx.fn}

\title{The product of x and log(x)}

\description{

  This function extends the product x.log(x) from positive real x to also
  include x = 0. 
  
}

\usage{

xlogx.fn(x)

}

\arguments{

  \item{x}{any scalar, vector or matrix x}
    
}

\details{

  The value of this function becomes zero at x = 0. It is used validly
  inside \code{clustglm} since the limit as x tends to zero of x.log(x)
  is zero. The function is mainly used internally in \code{clustglm}.
  
}

\value{

  A non-negtive real scalar, vector or matrix.
  
}

