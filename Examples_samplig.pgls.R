### Example: sampling.pgls
### Required packages:
library(caper);library(phylolm);library(phytools)
library(ggplot2)

set.seed(111)
### Simulating tree
tree<-pbtree(n=107)    
### Simulating response variable with phylogenetic signal
y <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))  
### Simulating explanatory variable
x <- y + rnorm(107,mean(y),1)     
plot(tree)

### Including Species names
sp <- tree$tip.label               
regre <- data.frame(sp,y,x)   

### Organizing comparative data for pgls:
c.data <- comparative.data(data=regre,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

### Fiting lm and pgls regressions:
mod.pgls <- pgls(y ~ x, c.data,"ML")
mod.lm <- lm(y ~ x,regre)
summary(mod.lm)
summary(mod.pgls)
  
### Regression Plots:
plot(y ~ x,pch=16)
abline(mod.lm,col="red",lwd=3)
abline(mod.pgls,col="blue",lwd=3)
legend("topleft",legend=c("LM","PGLS"),pch=c(16,16),col=c("red","blue"))
     
### Perform sensitive analysis:
result <- sampling.pgls(y~x,data=regre,phy=tree,times=2,breaks=seq(.1,.9,.1),lambda="ML",names.col="sp")
### You can specify the number of replicates per break interval:
result2 <- sampling.pgls(y~x,data=regre,phy=tree,times=50,breaks=c(.1,.5,.9),lambda="ML",names.col="sp")
