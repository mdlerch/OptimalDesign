\name{genCCD}
\alias{genCCD}
\title{Generate a central composite design}
\description{
    Generate a central composite design of arbitrary factors, with arbitrary
    number of center points, and user specified axial points.
}
\usage{genCCD(n_factors, n_center, alpha)}
\arguments{
    \item{n_factors}{an integer.  Number of factors.}
    \item{n_center}{an integer.  Number of center points.  Default 1}
    \item{alpha}{a number.  The radial distance from the center for the axial
        points.  Default, \code{sqrt(n_factors)}}
}
\details{
    Generate a central composite design of arbitrary number of factors with
    arbitrary number of center points, and user specified axial points.  The
    user can specify the radial distance of the axial points.
}
\value{A dataframe with \code{n_factors} columns.}
\author{Michael Lerch and John Sherrill}
\seealso{genFactorial}
\examples{
    genCCD(4)
    genCCD(4, 2)
    genCCD(4, 2, sqrt(2))
}
\keyword{design generation}
