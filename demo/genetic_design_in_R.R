geneticDesign <- function(formula, n, iterations = 1000, M = 10)
{

    parents <- list()
    for (i in 1:M)
    {
        parents[[paste(i)]] <- cbind(runif(n, -1 , 1), runif(n, -1, 1))
    }

    parent_eff <- numeric()
    for (i in 1:M)
    {
        parent_eff[i] <- getEff(formula, data.frame(parents[[paste(i)]]), criteria = "D")$D
    }

    alpha_blend <- runif(M, 0, 1)
    alpha_creep <- runif(M, 0, 1)
    alpha_mutat <- runif(M, 0, 1)
    alpha_bound <- runif(M, 0, 1)

    for (iter in 1:iterations)
    {
        second_parent <- sample(x = 1:M, size = M, replace = FALSE)

        children <- parents
        children <- geneticblend(children, parents, second_parent, alpha = alpha_blend)
        # cbind(children[["1"]], parents[["1"]], parents[[paste(second_parent[1])]])
        children <- geneticcreep(children, alpha = alpha_creep)
        # cbind(children[["1"]], parents[["1"]], parents[[paste(second_parent[1])]])
        children <- geneticmutat(children, alpha = alpha_mutat)
        # cbind(children[["1"]], parents[["1"]], parents[[paste(second_parent[1])]])
        children <- geneticbound(children, alpha = alpha_bound)
        cbind(children[["1"]], parents[["1"]], parents[[paste(second_parent[1])]])


        child_eff <- numeric()
        for (i in 1:M)
        {
            child_eff[i] <- getEff(formula, data.frame(children[[paste(i)]]), criteria = "D")$D
            if (!is.finite(child_eff[i]))
            {
                child_eff[i] <- 0
            }
        }

        print(cbind(parent_eff, child_eff))

        swap <- child_eff > parent_eff

        if (sum(swap) > 0)
        {
            mu_alpha_blend <- mean(alpha_blend[swap])
            mu_alpha_creep <- mean(alpha_creep[swap])
            mu_alpha_mutat <- mean(alpha_mutat[swap])
            mu_alpha_bound <- mean(alpha_bound[swap])

            alpha_blend[!swap] <- rnorm(length(sum(swap)), mu_alpha_blend, .1)
            alpha_creep[!swap] <- rnorm(length(sum(swap)), mu_alpha_creep, .1)
            alpha_mutat[!swap] <- rnorm(length(sum(swap)), mu_alpha_mutat, .1)
            alpha_bound[!swap] <- rnorm(length(sum(swap)), mu_alpha_bound, .1)

            alpha_blend[alpha_blend > 1] <- 1
            alpha_blend[alpha_blend < 0] <- 0
            alpha_creep[alpha_creep > 1] <- 1
            alpha_creep[alpha_creep < 0] <- 0
            alpha_mutat[alpha_mutat > 1] <- 1
            alpha_mutat[alpha_mutat < 0] <- 0
            alpha_bound[alpha_bound > 1] <- 1
            alpha_bound[alpha_bound < 0] <- 0
        }

        for (i in 1:M)
        {
            if (swap[i])
            {
                parents[[paste(i)]] <- children[[paste(i)]]
                parent_eff[i] <- child_eff[i]
            }
        }
    }

    parents
}

geneticblend <- function(children, parents, second_parent, alpha = NULL)
{
    if (is.null(alpha))
    {
        alpha <- rep(.2, length(children))
    }
    M <- length(children)
    nrows <- nrow(children[["1"]])
    for (child in 1:M)
    {
        for (r in 1:nrows)
        {
            if (runif(1) < alpha[child])
            {
                children[[paste(child)]][r, ] <- parents[[paste(second_parent[child])]][r, ]
            }
        }
    }
    children
}

geneticcreep <- function(children, size = .1, alpha = NULL)
{
    if (is.null(alpha))
    {
        alpha <- rep(.2, length(children))
    }
    M <- length(children)
    nrows <- nrow(children[["1"]])
    ncols <- ncol(children[["1"]])
    for (child in 1:M)
    {
        for (r in 1:nrows)
        {
            for (c in 1:ncols)
            {
                if (runif(1) < alpha[child])
                {
                    children[[paste(child)]][r, c] <- rnorm(1, 0, size) +
                        children[[paste(child)]][r, c]
                    children[[paste(child)]][which(children[[paste(child)]] >  1,
                                                   arr.ind = TRUE)] <- 1
                    children[[paste(child)]][which(children[[paste(child)]] < -1,
                                                   arr.ind = TRUE)] <- 0
                }
            }
        }
    }
    children
}

geneticmutat <- function(children, alpha = NULL)
{
    if (is.null(alpha))
    {
        alpha <- rep(.2, length(children))
    }
    M <- length(children)
    nrows <- nrow(children[["1"]])
    ncols <- ncol(children[["1"]])
    for (child in 1:M)
    {
        for (r in 1:nrows)
        {
            for (c in 1:ncols)
            {
                if (runif(1) < alpha[child])
                {
                    children[[paste(child)]][r, c] <- runif(1, -1, 1)
                }
            }
        }
    }
    children
}

geneticbound <- function(children, alpha = NULL)
{
    if (is.null(alpha))
    {
        alpha <- rep(.2, length(children))
    }
    M <- length(children)
    nrows <- nrow(children[["1"]])
    ncols <- ncol(children[["1"]])
    for (child in 1:M)
    {
        for (r in 1:nrows)
        {
            for (c in 1:ncols)
            {
                if (runif(1) < alpha[child])
                {
                    children[[paste(child)]][r, c] <- 1 * sign(children[[paste(child)]][r, c])
                }
            }
        }
    }
    children
}
