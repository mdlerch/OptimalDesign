library(OptimalDesign)

test <- function(ans, value, bigger = TRUE, context = "")
{
    m <- 1
    if (!bigger)
    {
        m <- -1
    }

    ans <- ans * m
    value <- value * m

    if (bigger)
    {
        ineq <- ">"
    } else
    {
        ineq <- "<"
    }
    msg <- paste(context, ":" , ineq, value, ". Observed:", ans)
    if (ans < value)
    {
        stop(msg)
    } else
    {
        cat(paste(msg, "\n"))
    }

}

# Make a dataset
X1 <- seq(-1, 1, .01)
X2 <- seq(-1, 1, .01)
X <- expand.grid(X1, X2)
names(X) <- c("X1", "X2")
formula <- ~X1 + X2

test(getEff(formula, optimalDesign(formula, X, 4, "D", 1000000), criteria = "D")$D,
     99, bigger = TRUE, context = "2^2 factorial Monte Carlo D")

test(getEff(formula, optimalDesign(formula, X, 4, "A", 1000000), criteria = "A")$A,
     73, bigger = TRUE, context = "2^2 factorial Monte Carlo A")

# test(getEff(optimalDesign(formula, X, 4, "I", 1000000), criteria = "I")$I,
#      3, bigger = FALSE, context = "2^2 factorial Monte Carlo I")

test(getEff(formula, optimalDesign(formula, X, 4, "G", 100000), evaluation = X, criteria = "G")$G,
     96, bigger = TRUE, context = "2^2 factorial Monte Carlo G")


context <- "2^2 factorial with genetic"

test(getEff(formula, geneticdesign(formula, X, 4, "D", 1000000, 10), criteria = "D")$D,
     99, bigger = TRUE, context = "2^2 factorial Monte Carlo D")

test(getEff(formula, geneticdesign(formula, X, A, "D", 1000000, 10), criteria = "A")$A,
     99, bigger = TRUE, context = "2^2 factorial Monte Carlo A")
