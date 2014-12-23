# From candidate, X, get intial design based on points
# with largest leverage
initialDesign <- function(X, size)
{
    # get the leverages of the design points (C++)
    l2norm <- initDesign(X)
    # order them (Move this to C++ at some point)
    order(l2norm, decreasing = TRUE)[1:size]
}
