\name{geneticdesign}
\alias{geneticdesign}
\title{Generate optimal designs using a genetic algorithm}
\description{
    Generate an optimal design on [-1,1] for an arbitrary linear model based on
    one of several design criteria.
}
\usage{geneticdesign(formula, dataframe, n, criterion = "D", iter = 100000, M = 10)}
\arguments{
    \item{formula}{an object of class 'formula'.  The model to be fitted based
        on the experiment}
    \item{dataframe}{a placeholder object of type 'dataframe'.  The dataframe
        must have a column for each variable in formula.  The contents of each
        column do not matter}
    \item{n}{an integer.  The number of design points to be used in the final
        design.}
    \item{criterion}{A string.  The criterion to maximize.  Currently, only "D"
        is supported}
    \item{iter}{an integer.  The number of iterations to use in the
        maximization routine.}
    \item{M}{an integer.  The number of parent or seed designs to be used in
        the algorithm}
}
\details{
    A genetic algorithm is used to produce a design with optimal
    characteristics.  All variables are allowed to range from -1 to 1.
}
\value{A 'data.frame'.}
\references{Use \url{} of \href{}{} for web}
\note{
   Not necessary, Pie charts are very bad way of ....
}
\author{Michael Lerch and John Sherrill}
\seealso{Point to code \code{\link{...}}}
\examples{
    geneticdesign(~X1 + X2, data.frame(X1=1, X2=1), "D")
}
\keyword{Optimal Design, experimental design, genetic algorithm}
