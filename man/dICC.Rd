
\name{dICC}
\alias{dICC}


\title{Distance-based Intra-Class Correlation Coefficient}

\description{Distance-based Intra-Class Correlation Coefficient (ICC) is an extension of the traditional univariate ICC to multivariate case, where the relationship between
the multivariate measurements is summarized in a distance matrix. It quantifies the ability of a measurement method in reproducing the inter-sample relationship. }

\usage{
 dICC(dist.mat, strata)
}


		
\arguments{

  \item{dist.mat}{a symmetrical distance matrix between all the replicates (technical and biological replicates).}
  \item{strata}{a factor with each level corresponding to a biological replicate.}

}



\value{

  Function \code{dICC} returns a list with the following component:
  \item{ICC}{the distance-based ICC value.}


}


\references{
Chen, J. & Zhang, X. 2022.  dICC: Distance-based Intraclass Correlation Coefficient for Metagenomic Reproducibility Studies.
submitted.

}

\author{Jun Chen and Xianyang Zhang}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{dICC.SE.asympt}}, \code{\link{dICC.SE.bt}}
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
dICC(dist.mat, strata)


}


\keyword{distance}
\keyword{multivariate}


