\name{genFactorial}
\alias{genFactorial}
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

    Factorial designs are some of the oldest experimental designs and
    have been in use since at least the 19th century. It is common to
    denote factorial designs by use of exponential notation. For
    example, a factorial design with 3 factors, each having 2 levels is
    denoted as a \eqn{2^3} factorial design and one having 7 factors,
    each having 4 levels is denotes as a \eqn{4^7} factorial design.
}
\value{A data.frame with \code{n_factors} columns.}
\author{Michael Lerch and John Sherrill}
\examples{
    genFactorial(3, 4)
}
\keyword{design generation}
