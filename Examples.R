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
#Simulating data and phylogeny:
N <- 100 # Number of species

####Example with a single tree#########
tree<-pbtree(n=N)
### Simulating response variable with phylogenetic signal and some standard error
y <- data.frame(resp=rTrait(n=5, tree, model=c("lambda"),parameters=list(lambda=.8)))
std <- function(x) sd(x)/sqrt(length(x)) #function for standard error
y$mean.resp<-rowMeans(y)
y$se.resp<-apply(y[,1:5],1,std)
### Simulating explanatory variable
x <- y$mean.resp + rnorm(N,mean(y$mean.resp),1)    
### Including Species names
sp <- tree$tip.label
#THE ANALYSIS
predtreeVar.pgls(resp=x,pred=y$mean.resp,se.pred=y$se.resp,tree,npred=2,method="normal",taxa.col=sp)


####Example with a set of trees of class multiphylo############
multitree <- rmtree(N=10,n=100)
### Simulating response variable with phylogenetic signal and some standard error
y <- data.frame(resp=rTrait(n=5, multitree[[1]], model=c("lambda"),parameters=list(lambda=.8)))
std <- function(x) sd(x)/sqrt(length(x)) #function for standard error
y$mean.resp<-rowMeans(y)
y$se.resp<-apply(y[,1:5],1,std)
### Simulating explanatory variable
x <- y$mean.resp + rnorm(N,mean(y$mean.resp),1)   
### Including Species names
sp <- multitree[[1]]$tip.label               
#THE ANALYSIS
predtreeVar.pgls(resp=x,pred=y$mean.resp,se.pred=y$se.resp,tree=multitree,npred=2,method="normal",taxa.col=sp)

