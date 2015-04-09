\name{genLatin}
\alias{genLatin}
\title{Generate a random latin square design}
\description{
    Generate a random latin square design of a particular size
}
\usage{genLatin(size, iterations)}
\arguments{
    \item{size}{an integer.  Size of latin square.}
    \item{iterations}{an integer. How long to run Markov chain.}
}
\details{
    A latin square can be defined as follows:

    Consider a square of dimension \eqn{n \times n} that within each cell,
    contains one of the first \eqn{n} latin letters. Say that each row
    contains one and only one appearance of each of these \eqn{n} letters
    and each column contains one and only one appearance of each of these
    \eqn{n} letters. This is a latin square.

    These squares may be used for experimental design when a small
    design is desired to study three factors that each have the same
    number of levels. Row, column, and letter are each considered a
    factor with \eqn{n} different levels.

    A single latin square of any size is essentially trivial to
    construct, however, generating \emph{random} latin squares is very
    difficult. There currently (April, 2015) is no closed form
    expression for the generation of latin squares, nor is there a
    closed form for the number of latin squares of a given size. As
    evidence of the difficulty of determining the number of latin
    squares of a given size, it is currently unknown how many latin
    squares exist of size 12 or greater.

    While there is no current method to generate exactly uniformly distributed
    latin squares, a paper published in 1996 by Mark Jacobson and Peter
    Matthews introduced an algorithm that creates a random walk through
    the set of latin squares and thus, by ergodicity, supplies a method
    for generating, essentially random latin squares.
    
    See Jacobson and Matthews, 1996.
}
\value{A matrix with \code{size} columns and rows.}
\author{Michael Lerch and John Sherrill}
\seealso{genFactorial, genCCD}
\examples{
    genLatin(4, 3)
}
\keyword{design generation}