### Example: sampling.pgls
### Required packages:
library(caper);library(phylolm);library(phytools)

### Before: copy and run functions: 
# sampling.pgls(); 
# influential.pgls() 
# plot.power.pgls()
# powerCtime()

N <- 50 # Number of species
### Simulating tree
tree<-pbtree(n=N)    
plot(tree)
### Simulating response variable with phylogenetic signal
Ly <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))  
### Simulating explanatory variable
Lx <- Ly + rnorm(N,mean(Ly),1)     

### Including Species names:
sp <- tree$tip.label               
regre <- data.frame(sp,Ly,Lx)   

### Organizing comparative data for pgls:
comp.data <- comparative.data(data=regre,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

### Linear regression (PGLS):
mod0 <- pgls(Ly ~Lx, data=comp.data,"ML")
summary(mod0)

### Estimating time before performing simulations:
powerCtime(Ly ~ Lx,data=comp.data,times=100,breaks=seq(.1,.9,.1))

### Example: sampling.pgls
samp1 <- sampling.pgls(Ly ~ Lx,data=comp.data,times=50)
### You can specify the number of replicates and break intervals:
samp1 <- sampling.pgls(Ly ~ Lx,data=comp.data,times=50,breaks=c(.1,.9,.1))

### You can specify the number of replicates per break interval:
samp2 <- sampling.pgls(Ly ~ Lx,data=comp.data,times=100,breaks=c(.1,.5,.9))

### Example: influence.pgls
influ1 <- influence.pgls(Ly ~ Lx,data=comp.data)
### Estimated parameters:
influ1$results
### Most influential species:
influ1[[4]]

### Visualizing Results:
plot.power.pgls(samp1,method="sampling")
plot.power.pgls(samp2,method="sampling")
plot.power.pgls(influ1,method="influence")




###########Trait variability analysis example#############
library(caper);library(phylolm);library(phytools)

####Example with a single tree#########
N <- 100 # Number of species
tree<-pbtree(n=N)
### Simulating response variable with phylogenetic signal
Ly <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))
### Simulating explanatory variable and its standard error
Lx <- data.frame(xmean=Ly + rnorm(N,mean(Ly),1),xse=rep(0.5,100))
### Including Species names
sp <- tree$tip.label
#THE ANALYSIS
predtreeVar.pgls(resp=Ly,pred=Lx$xmean,se.pred=Lx$xse,tree,npred=2,method="normal",taxa.col=sp)


####Example with a set of trees of class multiphylo############
N <- 100 # Number of species
multitree <- rmtree(N=50,n=100)
### Simulating response variable with phylogenetic signal
Ly <- rTrait(n=1, multitree[[1]], model=c("lambda"),parameters=list(lambda=.8))
### Simulating explanatory variable and its standard error
Lx <- data.frame(xmean=Ly + rnorm(N,mean(Ly),1),xse=rep(0.5,100))
### Including Species names
sp <- multitree[[1]]$tip.label               
#THE ANALYSIS
predtreeVar.pgls(resp=Ly,pred=Lx$xmean,se.pred=Lx$xse,tree=multitree,ntree=2,npred=2,method="normal",taxa.col=sp)
