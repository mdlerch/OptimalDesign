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
  Consider a design of size \eqn{N} with model matrix \eqn{X} and
  information matrix \eqn{M = X'X} having an average scaled prediction variance
  of \eqn{\sigma_{avg}} and maximum scaled predication variance of
  \eqn{\sigma_{max}} (see the vignette for the \code{OptimalDesign}
  package for a presentation on scaled prediction variance, design
  criterions, and other boring stuff). The criterion efficiencies are
  calculated with the following formulas:

  \deqn{\mbox{D-efficiency} = 100 \times \frac{|M|^{1/p}}{N}}
  \deqn{\mbox{A-efficiency} = 100 \times \frac{p/N}{\mbox{tr}(M^{-1})}}
  \deqn{\mbox{I-efficiency} = 100 \times \sigma_{avg}}
  \deqn{\mbox{G-efficiency} = 100 \times \sqrt{\frac{p/N}{\sigma_{max}}}.}
  
  The formulas used to compute the design efficiencies for the D, A, and
  G criterions are those provided by the SAS documentation for PROC
  OPTEX. The D- and A-efficiencies are the relative number of runs
  (expressed as percents) required by a hypothetical orthogonal design
  to achieve the same and , respectively; refer to Mitchell (1974b). "D"
  and "A" efficiencies are exact.  "I" and "G" efficiencies are 
  approximated based on the evaluation data frame.  Beware that approximate
  comparisons between "I" and "G" may be misleading.  See JOBO.
}
\value{A list of numbers}
\author{Michael Lerch and John Sherrill}
\examples{
     candidate <- genFactorial(100, 2)
     candidate <- candidate / max(candidate)
     getEff(~ X1 + X2, genFactorial(2, 2), candidate)
}
\keyword{zero or more keywords}
