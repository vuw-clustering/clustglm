## randpart.R

## Generates a random partition of n objects into c non-empty classes

randpart <- function(n,c){
    # n = number of objects
    # c = number of classes
    pp <- matrix(0,n,c)
    pp[1:c,] <- diag(rep(1,c))  # Ensures non-empty classes
    pp[(c+1):n,1] <- 1
    for (i in (c+1):n)
        pp[i,] <- sample(pp[i,])
    pp <- pp[sample(1:n),]
    pp
}
