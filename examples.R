### Example: power.pgls

### Required packages:
library(caper);library(phytools);library("phylolm")
set.seed(111)
### Simulating tree
tree<-pbtree(n=100)    
### Simulating response variable with phylogenetic signal
y <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))  
### Simulating explanatory variable
x <- varY + rnorm(100,mean(varY),1)     

### Including Species names
sp <- tree$tip.label               
data <- data.frame(sp,y,x,x2)   

### Organizing comparative data for pgls:
c.data <- comparative.data(data=data,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

### Fiting lm and pgls regressions:
mod.pgls <- pgls(y ~ x, c.data,"ML")
mod.lm <- lm(y ~ x,data)
summary(mod.lm)
summary(mod.pgls)
  
### Regression Plots:
plot(y ~ x,pch=16)
abline(mod.lm,col="red",lwd=3)
abline(mod.pgls,col="blue",lwd=3)
legend("topleft",legend=c("LM","PGLS"),pch=c(16,16),col=c("red","blue"))
     
### Perform sensitive analysis:
result <- power.pgls(y~x,data=data,phy=tree,times=10,breaks=seq(.1,.9,.1),lambda="ML",names.col="sp")
result2 <- power.pgls(y~x,data=data,phy=tree,times=50,breaks=c(.1,.5,.9),lambda="ML",names.col="sp")

### Visualizing results:
med <- with(result,tapply(betas,n.removs,mean))
Sdev <- with(result,tapply(betas,n.removs,sd))

n.re <- as.numeric(rownames(med))
n.sp <- nrow(c.data$data)-n.re 
result.med <- data.frame(med,Sdev,n.sp,n.re)
beta.mod.pgls <- summary(mod.pgls)[[c(5,2)]] 

ggplot(result.med,aes(y=med,x=n.sp))+
          geom_point(size=3,alpha=.7)+
          scale_x_continuous(breaks=seq(5, nrow(c.data$data), 5))+
          geom_errorbar(aes(ymin=med-Sdev, ymax=med+Sdev), width=.1)+
          geom_hline(yintercept= beta.mod.pgls,linetype=2,color="red")+
          xlab("Number of Species") + ylab("Estimated Beta")+
          geom_point(aes(x=nrow(c.data$data),y=beta.mod.pgls,size=3,colour="red"))+
          theme(legend.position="none",
                axis.text=element_text(size=14),       
                axis.title=element_text(size=20))

