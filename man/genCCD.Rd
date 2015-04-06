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

    Central Composite Designs were developed by Box and Wilson
    (1951) and are often used when a second order model is being
    considered. The design contains three components:
    \enumerate{
      \item a full factorial design
      \item a set of center points (points of the form \eqn{(0,0,\ldots,0)}
      \item a set of axial points, (points of the form
      \eqn{(x_1,x_2,\ldots,x_p)} where \eqn{x_i=\pm \alpha} for only one
      \eqn{i} and \eqn{x_j=0} all other \eqn{j} not equal to \eqn{i})
    }

    Let \eqn{F} denote the size of the factorial component of the
    design and let \eqn{T = 2k + n_{center}} the number of axial and
    center points in the design. To generate a rotatable design chose
    \eqn{\alpha = F^{1/4}} and to generate an orthogonal central
    composite design chose \eqn{\alpha = (F/4 (\sqrt{F+T} -
    \sqrt{F})^2)^{1/4}} (Myers, 2014).
  
    It's possible to use fractional factorial designs instead of full
    factorial designs in creating centeral composite designs. These are
    referred to as small composite designs (SCD) and where introduced by
    Hartley (1959). However, these are not yet supported by
    \code{OptimalDesign}.
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
