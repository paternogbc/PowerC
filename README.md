PowerC
======

PowerC is a a group of R functions to perform sampling effort analysis in comparative methods

###Functions:

	sampling.pgls - Simple sampling effort analysis for PGLS linear regression
	influential.pgls - Simple (leave-one-out deletion) diagnostics for PGLS linear regression
	plot.power.pgls - Plot results from sampling.pgls and influential.pgls



###Examples:
First copy and run these functions from the repository `PowerC`:   
	1. `sampling.pgls`  
          2. `influential.pgls`  
          3. `plot.power.pgls`  

**Required Packages:**
```{r}
library(caper);library(phytools);library("phylolm")
```

**Simulating data and phylogeny:**
```{r}
set.seed(111)
N <- 100 # Number of species
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
```

**Fitting regressions:**
```{r}
### Linear regression (PGLS):
mod0 <- pgls(y ~x, data=c.data,"ML")
summary(mod0)
```
**Regression Plots:**
```{r}
plot(y ~ x,c.data$data,pch=16)
abline(mod0,col="red",lwd=3)

```

**Performing Sensitive Analysis:** `sampling.pgls`
```{r}
samp1 <- sampling.pgls(y~x,data=regre,phy=tree,names.col="sp")
### You can specify the number of replicates per break interval:
samp2 <- sampling.pgls(y~x,data=regre,phy=tree,times=100,breaks=c(.1,.5,.9),names.col="sp")

```

**Performing influential analysis:** `influential.pgls`
```{r}
### Example: influence.pgls
influ1 <- influence.pgls(y ~ x,data=regre,phy=tree)
### Estimated parameters:
influ1$results
### Most influential species:
influ1[[4]]
```

**Visualizing Results:**
```{r,fig.show='hold'}
plot.power.pgls(samp1,method="sampling")
plot.power.pgls(samp2,method="sampling")
plot.power.pgls(influ1,method="influential")
```


