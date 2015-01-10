\name{getEff}
\alias{getEff}
\title{Get design efficiencies}
\description{
    Get design efficiencies for a specified design and model
}
\usage{getEff(formula, design, evaluation, criteria = c("D", "A", "I", "G"))}
\arguments{
    \item{formula}{A formula. Something like \code{~X1 + X2}}
    \item{design}{A data frame. The design points.  Columns should correspond
        to variables in the formula}
    \item{evaluation}{A data frame.  Optional. A set of points to evaluation
        the design criteria.  Necessary for "I" and "G" criteria.}
    \item{criteria}{A character or vector of characters. The design criteria to
        evaluate.}
}
\details{
    The specified design efficiencies are calculated.
}
\value{A list of numbers}
\note{
    "D" and "A" efficiencies are exact.  "I" and "G" efficiencies are
    approximated based on the evaluation data frame.  Be ware that approximate
    comparisons between "I" and "G" may be misleading.  See JOBO.
}
\author{Michael Lerch and John Sherrill}
\examples{
     candidate <- genFactorial(100, 2)
     candidate <- candidate / max(candidate)
     getEff(~ X1 + X2, genFactorial(2, 2), candidate)
}
\keyword{zero or more keywords}
