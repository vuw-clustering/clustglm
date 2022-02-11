## trymodels.R

## This file has code and comments for using "clustglm" to
## cluster rows and/or columns of a matrix of count data.

## Cut and paste sections of code into R.

rm(list=ls())

library("clustglm")

## Load the toy data set from P&A 2014

data(counts)

counts

class(counts)

dim(counts)


####################################################


## Set up the data frame:

?mat2df

counts.df <- mat2df(y = counts,
                    responsename = "Count",
                    factorname1 = "Row",
                    factorname2 = "Col")

str(counts.df)


## Fit some main-effects (ME) glms
## -------------------------------

## Use clustglm()rather than glm().
## Using glm() will give logLik, AIC and BIC out by a
## constant, so not comparable with later clustered models.
## clustglm() gives the true log likelihood, not adjusted by
## a constant, and so can be used to compare clustered and
## non-clustered models.

Null.out <- clustglm(formula = Count ~ 1,
                     data = counts.df,
                     family = "poisson")

RowME.out <- clustglm(formula = Count ~ Row,
                      data = counts.df,
                      family = "poisson")

ColME.out <- clustglm(formula = Count ~ Col,
                      data = counts.df,
                      family = "poisson")

RowPlusCol.out <- clustglm(formula = Count ~ Row + Col,
                           data = counts.df,
                           family = "poisson")

## Compare the models fitted so far:

?comparison

modellist <- list(Null = Null.out, Row = RowME.out,
                  Col = ColME.out, "Row+Col" = RowPlusCol.out)
comparison(modellist)

## Since the additive model (Row + Col) has the lowest AIC, we
## retain main effects of Row and Col in future clustered models.
## This is effectively standardising the data in order to focus
## on pattern detection.

## -----------------------------------------------------------
##                     Fit clustered models
## -----------------------------------------------------------

## Fit models with clustering of levels of factors Row and/or Col
## (i.e. row clustering, column clustering or biclustering of the
## original 8 by 10 matrix of counts).

## Here we retain the two main effects (Row and Col) from
## the unclustered models. This enables us to detect remaining
## patterns which occur in the interaction term.
## term, giving a broad view of patterns for ease of
## interpretation and plotting.



## --------------------------------------------------------------
##                    Cluster the rows
## --------------------------------------------------------------

## Two row clusters:
## -----------------

?clustglm

r2.out <- clustglm(formula = Count ~ Row + Col + RowClust:Col,
                       family = "poisson",
         	       data = counts.df,
         	       fact4clust = "Row", nclus = 2,
         	       clustfactnames = "RowClust",
         	       start.control = list(randstarts = 100),
         	       verbose = 1)

names(r2.out)
r2.out$LLint    ## -53.852
r2.out$AIC      ## 163.704

## Is this better than the unclustered models?
## Add this model to the model list, look for the lowest AIC.

modellist$r2 <- r2.out
comparison(modellist)


## Inspect posterior probabilities:
r2.out$pp.list   ## Posterior probabilities
round(r2.out$pp.list[[1]],3)

## Assign each A to its most likely cluster:
pp.mat <- r2.out$pp.list[[1]]
split(rownames(pp.mat), apply(pp.mat,1,which.max))

## $`1`
## [1] "A" "B" "C" "E" "G"
##
## $`2`
## [1] "D" "F" "H"

## There is a clear split between rows A,B,C,E,G and D,F,H.

## Profile plot
## ------------

## Profile plots and ordination plots use the interaction parameter
## estimates as plotting coordinates. See:
findpars(r2.out)
## which uses weighted sum-to-zero constraints.

## The following command does the profile plot and saves the coordinates

r2.profile <- profileplot(model = r2.out,
            x.factor = "Col",
            trace.factor = "RowClust",
	    main = "Two row clusters",
            legend = TRUE, xlab = "Col")

## Sorting the Col by gamma value can give a clearer picture:

r2.profile.sort <- profileplot(model = r2.out,
            x.factor = "Col",
            trace.factor = "RowClust",
	    main = "Count data: Two row clusters, columns ordered",
            sort.x = 2,  ## Sort for ascending values in cluster 2
            legend = TRUE, xlab = "Col")

## Now columns with similar patterns are adjacent.
## Interaction coefficients near zero are not showing patterns which are
## strongly marked. However, more extreme values suggest that column d
## is strongly associated with row cluster 1, while columns a, c, e and h
## are strongly associated with row cluster 2, with counts higher than
## would be expected from an additive (main effects) model.

## Ordination plot
## ---------------

## Two row clusters provide an ordination of the columns in 1D. The
## dimension reduction is caused by the weighted sum-to-zero parameter
## constraints in the model (linearity in the link scale).

r2.ord <- ordplot(model = r2.out,
        ord.factor = "Col",
        clust.factor = "RowClust",
	main = "Column ordination by two row clusters")

## The line has been swung around and displayed vertically. Both the
## sorted profile plot and the ordination plot show the same ordering
## of the columns, but the ordination plot also shows the spacing
## between them, with closer pairs being more alike in their patterns.

## For an ordination in two dimensions, fit a model with three row clusters.

## Three row clusters:
## -------------------

r3.out <- clustglm(formula = Count ~ Row + Col + RowClust:Col,
                       family = "poisson",
         	       data = counts.df,
         	       fact4clust = "Row", nclus = 3,
         	       clustfactnames = "RowClust",
         	       start.control = list(randstarts = 100),
         	       verbose = 0)

r3.out$LLint    ##  -40.419
r3.out$AIC      ##   158.839   An improvement over 2 row clusters

r3.out$pp.list
round(r3.out$pp.list[[1]],3)

pp.mat <- r3.out$pp.list[[1]]
## Find the estimated proportions, pi:
apply(pp.mat,2,mean)
## Find the effective size of each cluster:
apply(pp.mat,2,sum)
## Assign each row to its most likely cluster:
split(rownames(pp.mat), apply(pp.mat,1,which.max))

## $`1`
## [1] "B" "G"
##
## $`2`
## [1] "A" "C" "E"
##
## $`3`
## [1] "D" "F" "H"

## Profile plot:

r3.profile.sort <- profileplot(model = r3.out,
            x.factor = "Col",
            trace.factor = "RowClust",
	    main = "Three Row clusters, Col sorted using row cluster 3",
            sort.x = 3,  ## Sort Col for ascending values in cluster 3
            legend = TRUE, xlab = "Col")

## Clusters 2 and 3 each contain effectively three Rows, and so give more
## stable traces than cluster 1, with effectively only two rows.
## Cluster 1 is smaller and shows more fluctuation.

## Ordination plot:

## With three row clusters the 3 by J matrix of interaction parameters give
## a 2D ordination of the columns. See
findpars(r3.out)

## The ordination plot is

r3.ord <- ordplot(model = r3.out,
        ord.factor = "Col",
        clust.factor = "RowClust",
	main = "Column ordination by three row clusters")

## Each point represents a column. The points are in 3D, but are coplanar
## because of linearity in the link space. The plane has been rotated down
## to show the points in an x-y plane. The distances apart on this plot are
## consistent with the sorted profile plot and ordination plot from the
## model with two row clusters. This plot shows a "horseshoe effect" similar
## to those in matrix-based ordination.


## --------------------------------------------------------------
##                    Cluster the columns
## --------------------------------------------------------------

## Two column clusters:
## --------------------

c2.out <- clustglm(formula = Count ~ Row + Col + Row:ColClust,
                     family = "poisson",
                     data = counts.df,
                     fact4clust = "Col", nclus = 2,
                     clustfactnames = "ColClust",
                     start.control = list(randstarts = 100),
                     verbose = 1)

names(c2.out)
c2.out$LLint   ##  -52.157    A considerable improvement
c2.out$AIC     ##  156.314    over the main effects model

## See posterior probabilities:
c2.out$pp.list   ## Posterior probabilities
round(c2.out$pp.list[[1]],3)

## Assign each B to its most likely cluster:
pp.mat <- c2.out$pp.list[[1]]
split(rownames(pp.mat), apply(pp.mat,1,which.max))

## $`1`
## [1] "a" "c" "e" "h"
##
## $`2`
## [1] "b" "d" "f" "g" "i" "j"

## The clusters may come out the other way around. Labels can switch.

## Do a profile plot (interaction plot) of the interaction term, having
## allowed for the main effects. Sort rows to show increasing interaction
## coefficients for the larger column cluster.

c2.profile.sort <- profileplot(model = c2.out,
            x.factor = "Row",
            trace.factor = "ColClust",
	    main = "Profile plot: Two column clusters",
	    sort = 2,
            legend = TRUE, xlab = "Row")

## Do an ordination in 1D:

c2.ord <- ordplot(model = c2.out,
        ord.factor = "Row",
        clust.factor = "ColClust",
	main = "Row ordination by two column clusters")

## Three column clusters:
## -----------------------

c3.out <- clustglm(formula = Count ~ Row + Col + Row:ColClust,
                       family = "poisson",
         	       data = counts.df,
         	       fact4clust = "Col", nclus = 3,
         	       clustfactnames = "ColClust",
         	       start.control = list(randstarts = 100),
         	       verbose = 1)

names(c3.out)
c3.out$LLint   ##  -41.355
c3.out$AIC     ##  152.71

## See posterior probabilities:
c3.out$pp.list   ## Posterior probabilities
round(c3.out$pp.list[[1]],3)

## Assign each B to its most likely cluster:
pp.mat <- c3.out$pp.list[[1]]
split(rownames(pp.mat), apply(pp.mat,1,which.max))

## $`1`
## [1] "a" "c" "e" "h"
##
## $`2`
## [1] "d" "f"
##
## $`3`
## [1] "b" "g" "i" "j"

## Do a profile plot of the three clusters:

c3.profile <- profileplot(model = c3.out,
            x.factor = "Row",
            trace.factor = "ColClust",
	    main = "Profile plot: Three Col clusters",
      	    legend = TRUE, xlab = "Row")

## See also:
findpars(c3.out)
apply(pp.mat,2,sum)

## A value of -12.5 dominates this plot. The small cluster c("d","f") has
## all zero counts in row H, so that with a log link the interaction
## term is trying to be -infinity. This can happen with small data sets
## or with sparse data sets with many zeros. A possible adjustment would
## replace the long line with a short line and arrow to indicate it is
## trying to go to -infinity. The rest of the plot expands to show the
## other patterns more clearly.

## Do an ordination plot, in which the rows are ordinated by the three
## column clusters.

c3.ord <- ordplot(model = c3.out,
        ord.factor = "Row",
        clust.factor = "ColClust",
	main = "Row ordination by three column clusters")

## This plot also shows the difficult caused by zeros for all of one column
## cluster in row H.

## ----------------------------------------------------------------------------
##              Model comparisons so far
## ----------------------------------------------------------------------------

mod.list <- list(Null = Null.out,
                 Row = RowME.out,
                 Col = ColME.out,
		 "Row+Col" = RowPlusCol.out,
		 r2 = r2.out,
		 r3 = r3.out,
		 c2 = c2.out,
		 c3 = c3.out)

## Do the comparison:

mod.table <- comparison(mod.list)

##             logL npar     AIC relAIC     BIC relBIC
## Null    -105.375    1 212.750 60.040 215.132  0.000
## Row      -93.486    8 202.972 50.262 222.028  6.896
## Col      -94.934   10 209.867 57.157 233.687 18.555
## Row+Col  -83.044   17 200.089 47.379 240.583 25.451
## r2       -53.852   28 163.704 10.994 230.401 15.269
## r3       -40.420   39 158.839  6.129 251.738 36.606
## c2       -52.157   26 156.314  3.604 218.247  3.115
## c3       -41.355   35 152.710  0.000 236.081 20.949

## "relAIC" is the value relative to the minimum.
## The model with lowest AIC has relAIC = 0.  Similarly BIC.
## AIC selects the model with three column clusters, while BIC, with its
## higher penalty on the number of parameters, selects th null model.



## ----------------------------------------------------------------------
##       Biclustering: Simultaneous clustering of rows and columns
## ----------------------------------------------------------------------

## Two row clusters and column clusters
## ------------------------------------

r2c2.out <- clustglm(formula = Count ~ Row + Col + RowClust:ColClust,
                      family = "poisson",
                      data = counts.df,
                      fact4clust = c("Row","Col"), nclus = c(2,2),
                      clustfactnames = c("RowClust","ColClust"),
                      start.control = list(randstarts = 100),
                      verbose = 1)

r2c2.out$LLint    ## -69.359
## Tried again, better answer.
r2c2.out2$LLint   ## -66.905
## Third try, same as second.
r2c2.out3$LLint   ## -66.905

## Try allocated start using single-mode clustering:
pp.start <- list(r2.out$pp.list[[1]],c2.out$pp.list[[1]])

r2c2.out4 <- clustglm(formula = Count ~ Row + Col + RowClust:ColClust,
                      family = "poisson",
                      data = counts.df,
                      fact4clust = c("Row","Col"), nclus = c(2,2),
                      clustfactnames = c("RowClust","ColClust"),
                      start.control = list(alloc = pp.start),
                      verbose = 1)
r2c2.out4$LLint   ## -66.905, same as above.

## The biclustered model has fewer parameters than the single-mode
## clustering, and cannot achieve such a high log likelihood.
## Update using this fit:

r2c2.out <- r2c2.out4

## Three row clusters and column clusters
## --------------------------------------

pp.start <- list(r3.out$pp.list[[1]],c3.out$pp.list[[1]])

r3c3.out <- clustglm(formula = Count ~ Row + Col + RowClust:ColClust,
                      family = "poisson",
                      data = counts.df,
                      fact4clust = c("Row","Col"), nclus = c(3,3),
                      clustfactnames = c("RowClust","ColClust"),
                      start.control = list(alloc = pp.start),
                      verbose = 1)

r3c3.out$LLint   ## -60.713, higher than r2c2

pp1 <- round(r3c3.out$pp.list[[1]],3)
pp2 <- round(r3c3.out$pp.list[[2]],3)
pp1
pp2
## Allocate to most likely clusters:

split(rownames(pp1), apply(pp1,1,which.max))

## $`1`
## [1] "B" "G"
##
## $`2`
## [1] "A" "C" "E"
##
## $`3`
## [1] "D" "F" "H"

split(rownames(pp2), apply(pp2,1,which.max))

## $`1`
## [1] "a" "c" "e" "h"
##
## $`2`
## [1] "d" "f" "g" "i"
##
## $`3`
## [1] "b" "j"



## ---------------------------------------------------------
##      Find parameters with natural constraints
## ---------------------------------------------------------

r3c3.pars <- findpars(r3c3.out)
r3c3.pars

## $`Overall wt.mean`
## [1] 1.688166
##
## $Row
##           A           B           C           D           E
##  0.02683895  0.26980105 -0.31266806 -0.22467082 -0.08528024
##            F           G           H
## 	     0.35338368  0.17573975 -0.20314431
##
## $Col
##           a           b           c           d           e
## -0.22094867  0.13529454 -0.02033239 -0.05251558 -0.38584907
##           f           g           h           i           j
##  0.03264256  0.16617244 -0.24864924  0.21771495  0.37647046
##
## $`RowClust:ColClust`
##            ColClust1  ColClust2  ColClust3
## RowClust1 -0.1486823 -0.1027084  0.4698188
## RowClust2 -0.4831763  0.5472480 -0.1225665
## RowClust3  0.5812109 -0.4786263 -0.1888967

## The parameter estimates above are from the fitted values on the link scale;
## since this glm family is Poisson, they are on the log scale.
## Each parameter is a weighted mean of fitted values, the weight being the
## expected number of observations contributing to its estimation.

## The parameters are in a list with four components:
## (i) "nu", the intercept, is the overall weighted mean fitted value,
## (ii) "alpha" is a vector of deviations from "nu", one for each row of
##      the original matrix,
## (iii) "beta" is a vector of deviations  from "nu", one for each column of
##      the original matrix, and
## (iv) "gamma" is the interaction parameter; for each combination of the
##      interaction factor levels it is the amount of fitted value remaining
##      after allowing for the intercept and main effects.

## For this balanced data set the alphas sum to zero, since each level of the Row
## factor has the same number of observations, so the weights are equal; similarly
## the betas add to zero. We check this for this example:

sum(r3c3.pars[[2]])   ## 4.739e-10
sum(r3c3.pars[[3]])   ## 2.748e-10

## The row and column constraints on the gamma matrix are for weighted sums to
## be zero, the expected number of observations in row cluster r and column
## cluster c being N.pi(r).kappa(c). Factoring out N (number of observations),
## we check that the matrix with terms pi(r).kappa(c).gamma(r,c) has row and
## column sums zero:

pi.vect <- r3c3.out$pi.list[[1]]
kappa.vect <- r3c3.out$pi.list[[2]]
wted.gamma.mat <- diag(pi.vect) %*% r3c3.pars[[4]] %*% diag(kappa.vect)
round(wted.gamma.mat,3)
apply(wted.gamma.mat,1,sum)
apply(wted.gamma.mat,2,sum)

## These values are zero to within rounding errors of three or four decimal places.

## To interpret the patterns after allowing for main efects, look at the gamma
## estimates:

gam33.mat <- round(r3c3.pars[[4]],2)
gam33.mat

##           ColClust1 ColClust2 ColClust3
## RowClust1     -0.15     -0.10      0.47
## RowClust2     -0.48      0.55     -0.12
## RowClust3      0.58     -0.48     -0.19

## Strongest positive associations are between:
## (i) RowClust3 and ColClust1, i.e. (D,F,H) and (a,c,e,h)
## (ii) RowClust2 and ColClust2, i.e. (A,C,E) and (d,f,g,i)
## (iii) RowClust1 and ColClust3, i.e. (B,G) and (b,j)

## Strongest negative associations are between:
## (i) RowClust2 and ColClust1, i.e. (A,C,E) and (a,c,e,h)
## (ii) RowClust3 and ColClust2, i.e. (D,F,H) and (d,f,g,i)


## --------------------------------------------------
##        Update the model comparison table:
## --------------------------------------------------

mod.list$r2c2 <- r2c2.out
mod.list$r3c3 <- r3c3.out

## Do the comparison:

mod.table <- comparison(mod.list)

##             logL npar     AIC relAIC     BIC relBIC
## Null    -105.375    1 212.750 60.040 215.132  0.000
## Row      -93.486    8 202.972 50.262 222.028  6.896
## Col      -94.934   10 209.867 57.157 233.687 18.555
## Row+Col  -83.044   17 200.089 47.379 240.583 25.451
## r2       -53.852   28 163.704 10.994 230.401 15.269
## r3       -40.420   39 158.839  6.129 251.738 36.606
## c2       -52.157   26 156.314  3.604 218.247  3.115
## c3       -41.355   35 152.710  0.000 236.081 20.949
## r2c2     -66.905   22 177.810 25.100 230.215 15.083
## r3c3     -60.713   29 179.426 26.716 248.505 33.373

## Note: mod.table is actually a data frme, see
class(mod.table)



##############################################################
##                 Biplot from biclustering
##############################################################

## The function "clustbiplot" is not yet fully implemented

