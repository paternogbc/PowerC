\name{influence.pgls}
\alias{influence.pgls}

\title{
Influence analysis for PGLS 
}

\description{
This function performs 'single deletion analysis' for pgls linear regression
}
\usage{
influence.pgls(formula, data, lambda = "ML")
}
\arguments{
  \item{formula}{
A model formula
}
  \item{data}{
A 'comparative.data' object containing the covariance matrix and data to be used in the model (see function comparative.data in the caper package)}
  \item{lambda}{
A value for the lambda transformation. If NULL, lambda="ML"}
}
\details{
This function can only deal with simple linear regression y = a bx + e. Future implementation will deal with more complex models. 
}
\value{
Output is a list with the following components:
\item{errors}{Character vector with species that showed error during pgls fit}
\item{formula}{The model formula}
\item{original_model_estimates}{Original estimates from the full model}
\item{original_beta_95_IC}{Beta confidence interval from original model }
\item{most_influential_species}{The five most influential species for beta and intercept estimation}
\item{results}{A data frame with: Species name, estimated beta, difference from original beta (DFbeta), estimated intercept, difference from original intercept (DFintercept) and p-values for each simulation (single deletion)}
}

\references{
R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic analysis and comparative data: A test and review of evidence. American Naturalist, 160:712-726, 2002.

}
\author{
Gustavo Brant Paterno; Gijsbert Werner }
\note{
}

\section{Warning }{This code is not fully checked, please be aware} 

\seealso{
\code{\link{pgls}}
}
\examples{
### Required packages:
library(caper);library(phylolm);library(phytools)

N <- 80 
tree<-pbtree(n=N)    
Ly <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))  
Lx <- Ly + rnorm(N,mean(Ly),1)     
sp <- tree$tip.label               
regre <- data.frame(sp,Ly,Lx)   
comp.data <- comparative.data(data=regre,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

### Performing influence analysis:
influ <- influence.pgls(Ly ~ Lx,data=comp.data)
### Estimated parameters:
influ$results
### Most influential species:
influ[[4]]

### Visualizing Results:
plot.power.pgls(influ1,method="influence")