summary.clustglm <- function (object, ...) {
  ans <- list(glmsummary = NextMethod("summary", object), pp.list = object$pp.list)  ## interim
  class(ans) <- "summary.clustglm"
  ans
}

print.summary.clustglm <- function (x, ...) {
  print(x$glmsummary)  ## cat('not ready yet\n')
  cat("\n")
  if (!is.null(x$pp.list)) {
      cat("Posterior probability of cluster membership\n")
      cat("-------------------------------------------\n")
      print(x$pp.list)  
  }
  invisible(x)
}

## just playing with options 2016-12-14
print.clustglm <- function(x, ...) {
  NextMethod("print",x)
  cat("\n")
  if (!is.null(x$pp.list)) {
      cat("Posterior probability of cluster membership\n")
      cat("-------------------------------------------\n")
      print(x$pp.list)  
  }
  invisible(x)
}
