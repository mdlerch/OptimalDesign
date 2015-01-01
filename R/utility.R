#
# Collection of small utility functions
#

getEff <- function(formula, design, candidate)
{
	candidate <- model.matrix(formula, candidate)
    X <- model.matrix(formula, design)
	p <- ncol(X)
	N <- nrow(X)
		
    # information matrix
    infoMat <- solve(t(X)%*%X)
	
    # evaluate efficiencies
    dEff <- 100 / (det(infoMat)^(1/p) * N)
    aEff <- 100*p / (sum(diag(infoMat)) * N)

	efficiencies <- list(dEff=dEff, aEff=aEff)
	
	# calculate iff candidate set provided
    # I don't know about these to be honest
	if(!missing(candidate)){
		spv <- apply(N * t(candidate) * (infoMat%*%t(candidate)), 2, sum)
		cat(spv)
	 	gEff <- 100 * sqrt(p/N) / max(spv^.5)
		iEff <- mean(spv^.5)
		
		efficiencies <- c(efficiencies, list(gEff=gEff, iEff=iEff))
    }
    
    efficiencies
}

#
# Generate factorial designs. It's recursive!
#
# rename these functions as you see fit (or the variables for that matter)
#

genFact <- function(n_levels, n_terms, i=n_terms){
    # center levels to zero and get scaling like in AlgDesign
	levels <- scale(1:n_levels - n_levels/2 -.5, scale=(n_levels%%2/2+.5))
	this <- rep(rep(levels ,n_levels^(n_terms-i)), rep(n_levels^(i-1), n_levels^(n_terms-i+1)))

	if(i>1) data.frame(unname(cbind( this, genFact(n_levels,n_terms,i-1) )))
	else this
}

genBBD <- function()
{
	
}

genCCD <- function()
{
	
}
