\name{spvEquation}
\alias{spvEquation}
\title{Print human readable SPV}
\description{
    Print human readale scaled prediction variance (SPV) equation
}
\usage{spvEquation(formula, design)}
\arguments{
    \item{formula}{A formula. Something like \code{~X1 + X2}}
    \item{design}{A data frame. The design points.  Columns should correspond
        to variables in the formula}
}
\details{
  Generates a human readable scaled prediction variance (SPV) function given a
  specified design \eqn{X} and formula. See the vignette for the
  \code{OptimalDesign} package for a presentation on scaled prediction
  variance function and it's utility.

  This function utilizes a homebrewed computer algebra system and the
  fundamental theorem of arithmetic to generate the SPV function.
}

\value{A string containing the SPV equation}
\note{It is recommended that the variable names X1, X2, X3, ... be used
    (or something similar). Long variable names will make the SPV
    equation hard to read.

    At the moment, this function assumes there is an intercept in the model.
}
\author{Michael Lerch and John Sherrill}
\examples{
     design <- genFactorial(10, 2)
     spvEquation(~., design)
}
\keyword{scaled prediction variance, experimental design}
