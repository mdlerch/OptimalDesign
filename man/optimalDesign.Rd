\name{optimalDesign}
\alias{optimalDesign}
\title{Generate optimal designs}
\description{
    Generate an optimal design from a candidate set of design points based on
    one of several design criteria.
}
\usage{optimalDesign(formula, candidate, n, criterion = "D", iter = 10000)}
\arguments{
    \item{formula}{an object of class 'formula'.  The model to be fitted based
        on the experiment.}
    \item{candidate}{a data frame. The candidate set of design points}
    \item{n}{an integer.  The number of design points to be used in the final
        design.}
    \item{criterion}{A string.  The criterion to maximize.  One of "D", "A",
        "G", or "I".}
    \item{iter}{an integer.  The number of iterations to use in the
        maximization routine.}
}
\details{
    For the "D" criterion, the Fedorov algorithm is used to find an optimal
    design.  The "A", "G", and "I" criteria are not supported yet.
}
\value{A data.frame.  A subset of \code{candidate} of length \code{n}.  Points
    may be repeated.}
\references{Use \url{} or \href{}{} for web}
\author{Michael Lerch and John Sherrill}
\examples{
    cand <- genFactorial(2, 4)
    optimalDesign(~ X1 + X2, cand, 4, "D")
}
