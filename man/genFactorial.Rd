\name{genFactorial}
\alias{genFactorial}
\title{Generate a factorial design}
\description{
    Generate a factorial design at arbitrary levels and variables
}
\usage{genFactorial(n_factors, n_levels)}
\arguments{
  \item{n_factors}{an integer.  The number of factors.}
  \item{n_levels}{an integer. The number of levels for each factor.}
}
\details{
    Generate a factorial design of arbitrary levels and factors.  The factor
    levels are scaled to be evenly spaced integers centered at 0.

    Factorial designs are some of the oldest experimental designs and
    have been in use since at least the 19th century. A factorial design
    is essentially a design in which all combinations of levels for each
    factor are included. For example, say an experiment has 3 factors
    (\eqn{X1}, \eqn{X2}, \eqn{X3}) and each factor has 2 levels, then
    the factorial design would contain the points
    \tabular{ccc}{
      X1 \tab X2 \tab X3 \cr
      -1 \tab -1 \tab -1 \cr
      -1 \tab -1 \tab 1 \cr
      -1 \tab 1 \tab -1 \cr
      -1 \tab 1 \tab 1 \cr
      1 \tab -1 \tab -1 \cr
      1 \tab -1 \tab 1 \cr
      1 \tab 1 \tab -1 \cr
      1 \tab 1 \tab 1
    }.
    It is common to denote factorial designs by use of exponential
    notation. For example, a factorial design with 3 factors, each
    having 2 levels, e.g. the above design, is denoted as a \eqn{2^3}
    factorial design and one having 7 factors, each having 4 levels is
    denotes as a \eqn{4^7} factorial design. See (Montgomery, 183) for more.

    It's common that full factorial designs are too large when
    constraints (logistic, monetary, etc.) are taken into
    consideration. In this case, often a \emph{fractional} factorial
    design can be created. These designs are not currently supported by
    \code{OptimalDesign} but are planned for future development.
}
\value{A data.frame with \code{n_factors} columns.}
\author{Michael Lerch and John Sherrill}
\examples{
    genFactorial(3, 4)
}
\keyword{design generation}
