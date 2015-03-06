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
    The scaled prediction variance equation is an expression for the
    variance of predicted values determined from a model fit using a
    particular design.

    You might want to check out
    \link{http://en.wikipedia.org/wiki/Mean_and_predicted_response}

    Probably should put in a nice discussion on prediction variance and
    the criteria generated.
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
