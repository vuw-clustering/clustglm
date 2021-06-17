# clustglm
Model-based clustering of binary, count and binomial data with pattern detection and covariates.

This package uses finite mixture models to perform single-mode clustering and bi/co-clustering of binary, count and binomial data. If you think of the data as a matrix with observations in the rows, and variables in the columns, single mode clustering is when you cluster just the observations (rows), or when you cluster just the variables (columns). Bi-clustering (also called co-clustering) is when you cluster both observations and variables simultaneously.

In an ecological dataset, the data might be presence/absence data, and the rows might be different species, and the columns might be different sites, and so you could cluster both species and sites to find out which groups of species occur at the same subset of sites.

This package incorporates the R function glm(), which allows you to use distribution families such as Bernoulli, Poisson and Binomial. You call the main function, `clustglm`, using the usual formula syntax, using additional arguments to specify whether you want to cluster the rows or the columns, but you can also include other covariates that are not part of the original data matrix. You can also specify a model in which there is a term for individual rows, as well as a term for the row clusters; i.e. unlike other clustering models, individual observations that are in the same cluster can have distinct patterns (species within a cluster can still have their own individual features).

Shirley Pledger (at Victoria University of Wellington, New Zealand) is the creator and author of this package. Co-authors: Murray Efford, Richard Arnold, Daniel Fernández, Ivy Liu, Lloyd Pledger. Louise McMillan is the maintainer.

If you need any help, please email louise.mcmillan@vuw.ac.nz.


## Citations
Please cite

Pledger, S., & Arnold, R. (2014). Multivariate methods using mixtures: Correspondence analysis, scaling and pattern-detection. Computational Statistics & Data Analysis, 71, 241-261.

````markdown
@article{pledger2014multivariate,
  title={Multivariate methods using mixtures: Correspondence analysis, scaling and pattern-detection},
  author={Pledger, Shirley and Arnold, Richard},
  journal={Computational Statistics \& Data Analysis},
  volume={71},
  pages={241--261},
  year={2014},
  publisher={Elsevier}
}
````
