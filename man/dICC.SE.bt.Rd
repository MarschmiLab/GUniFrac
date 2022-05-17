
\name{dICC.SE.bt}
\alias{dICC.SE.bt}


\title{Bootstrap Standard Error of Distance-based Intra-Class Correlation Coefficient}

\description{Calculate the standard error of the distance-based intra-class correlation coefficient based on bootstrapping. Biological replicates are sampled by replacement. 
More conservative than the asymptotic approach.}

\usage{
 dICC.SE.bt(dist.mat, strata, B = 199)
}


		
\arguments{

  \item{dist.mat}{a symmetrical distance matrix between all the replicates (technical and biological replicates).}
  \item{strata}{a factor with each level corresponding to a biological replicate. Supports an unequal number of technical replicates for each biological replicate.}
  \item{B}{integer, the number of bootstrap iterations.}
  

}


\value{

  Function \code{dICC} returns a list with the following component:
  \item{ICC}{the distance-based ICC value.}
   \item{SE}{the standard error.}


}


\references{
Chen, J. & Zhang, X. 2022.  dICC: Distance-based Intraclass Correlation Coefficient for Metagenomic Reproducibility Studies.
submitted.

}

\author{Jun Chen and Xianyang Zhang}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{dICC}}, \code{\link{dICC.SE.asympt}}
}

\examples{


# Generate the error-free measurements of 20 biological replicates, each with four dimensions
y <- matrix(rnorm(80), nrow = 20, ncol = 4)

# Generate two technical replicates for each biological replicate and add measurement errors
y1 <- y + matrix(rnorm(80, sd = 0.5), nrow = 20, ncol = 4)
y2 <- y + matrix(rnorm(80, sd = 0.5), nrow = 20, ncol = 4)
y12 <- rbind(y1, y2)

# Create the design vector
strata <- factor(rep(1 : 20, 2))

# Calculate the distance and distance-based ICC
dist.mat <- as.matrix(dist(y12))
dICC.SE.bt(dist.mat, strata)


}





\keyword{distance}
\keyword{multivariate}


