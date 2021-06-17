\name{findpars}
\alias{findpars}

\title{Find parameter estimates from the linear predictor}

\description{

  For a given model the estimated coefficients are reparameterized to
  show main effects and (optionally) interactions with sum-to-zero
  constraints. The interaction parameter, gamma, shows patterns after
  main effects are allowed for.
  
}

\usage{

   findpars(model)

}

\arguments{

  \item{model}{a model for which parameters with sum-to-zero constraints
    are required}
    
}

\details{

  For a model with linear predictor of the type y ~ A + B + C:D, the
  coefficients are reparameterized to have sum-to-zero constraints for
  unclustered factors, and weighted sum-to-zero constraints for clustered
  factors.
    
}

\value{

  A list with items nu = overall mean, alpha = a vector of deviations
  from nu for data with factor A at level a, beta = a vector of
  deviations from nu for data with factor B at level b, and gamma =
  deviation from the main-effects model nu + alpha(a) + beta(b). The
  gamma parameters give pattern detection information.

}

\seealso{

\code{\link{clustglm}}

}

\examples{

## To be provided.

}

\keyword{ parameter estimates }

