logLik.clustglm <- function (object, ...){
    if (is.null(object$fact4clust))
        ans <- object$LL
    if (is.character(object$fact4clust))
        ans <- object$LLint
    class(ans) <- "logLik.clustglm"
    ans
}


