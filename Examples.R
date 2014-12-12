### Example: sampling.pgls
### Required packages:
library(caper);library(phylolm);library(phytools)

### Before: copy and run functions: sampling.pgls(); influential.pgls() and plot.power.pgls()

set.seed(111)
N <- 50 # Number of species
### Simulating tree
tree<-pbtree(n=N)    
### Simulating response variable with phylogenetic signal
y <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))  
### Simulating explanatory variable
x <- y + rnorm(N,mean(y),1)     

### Including Species names
sp <- tree$tip.label               
regre <- data.frame(sp,y,x)   

### Organizing comparative data for pgls:
c.data <- comparative.data(data=regre,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

### Linear regression (PGLS):
mod0 <- pgls(y ~x, data=c.data,"ML")
summary(mod0)

### Example: sampling.pgls
samp1 <- sampling.pgls(y~x,data=regre,phy=tree,times=10,breaks=seq(.1,.9,.1),lambda="ML",names.col="sp")

### You can specify the number of replicates per break interval:
samp2 <- sampling.pgls(y~x,data=regre,phy=tree,times=50,breaks=c(.1,.5,.9),lambda="ML",names.col="sp")

### Example: influence.pgls
influ1 <- influence.pgls(y ~ x,data=regre,phy=tree,lambda="ML")
### Estimates values:
influ1$estimates
### Most influential species:
influ1[[5]]

### Visualizing Results:
plot.power.pgls(samp1,method="sampling")
plot.power.pgls(samp2,method="sampling")
plot.power.pgls(influ1,method="influential")
