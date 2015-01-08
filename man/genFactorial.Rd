\name{genFactorial}
\title{Generate a factorial design}
\description{
    Generate a factorial design at arbitrary levels and variables
}
\usage{genFactorial(n_levels, n_terms)}
\arguments{
    \item{n_levels}{an integer. The number of levels for each factor.}
    \item{n_factors}{an integer.  The number of factors.}
}
\details{
    Generate a factorial design of arbitrary levels and factors.  The factor
    levels are scaled to be evenly spaced integers centered at 0.
}
\value{A data.frame with \code{n_factors} columns.}
\author{Michael Lerch and John Sherrill}
\examples{
    \dontrun{genFactorial(3, 4)}
}
\keyword{design generation}
