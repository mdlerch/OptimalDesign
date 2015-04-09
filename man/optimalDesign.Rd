\name{optimalDesign}
\alias{optimalDesign}
\title{Generate optimal designs}
\description{
    Generate an optimal design from a candidate set of design points based on
    one of several design criteria.
}
\usage{optimalDesign(formula, candidate, n, criterion = "D", iter = 10000,
                          evaluation = candidate, repeated = TRUE)}
\arguments{
    \item{formula}{an object of class 'formula'.  The model to be fitted based
        on the experiment.}
    \item{candidate}{a data frame. The candidate set of design points}
    \item{n}{an integer.  The number of design points to be used in the final
        design.}
    \item{criterion}{a string.  The criterion to maximize.  One of "D", "A",
        "G", or "I".}
    \item{iter}{an integer.  The number of iterations to use in the
        maximization routine.}
    \item{evaluation}{a data frame.  The set of points on which to evaluate "G"
        or "I" criteria.  Defaults to candidate set}
    \item{repeated}{boolean.  If TRUE allow each candidate point to potentially
        repeated.  If FALSE, use each candidate point at most once}
}
\details{
    A Monte-Carlo-modified Fedorov search algorithm is used to find an optimal
    design with respect to a particular design criterion. Currently, the
    D, A, I, and G criterions are supported, but the G criterion does not
    work well. See the vignette for the \code{OptimalDesign} package for
    a presentation of these criterions and an explination of the
    modified Fedorov search algorithm.
}
\value{A data.frame.  A subset of \code{candidate} of length \code{n}.  Points
    may be repeated.}
\references{Use \url{} or \href{}{} for web}
\author{Michael Lerch and John Sherrill}
\examples{
    cand <- genFactorial(2, 4)
    optimalDesign(~ X1 + X2, cand, 4, "D")
}
