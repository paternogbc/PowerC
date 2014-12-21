#Gijsbert Werner, Vrije Universiteit Amsterdam
#Email: g.d.a.werner@vu.nl

#December 2014

#Create a R-function to analyse sensitivity to data uncertainty in binary trait evolution (specifically in corHMM)

library(phytools);library(phylolm);library(corHMM)

#Simulate a tree
set.seed(111)
N <- 100 # Number of species
### Simulating tree
tree<-pbtree(n=N)    
plot(tree)
### Simulating response variable with phylogenetic signal
Ly<-rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8)) #Use exact same simulations as previously as basis

#Everything above this, doesn't need to be redone, it's already in Gustavo's code
X<-cbind(rep(1,100),Ly)
y<-rbinTrait(n=1, phy=tree,beta=c(-1,0.5), alpha=1, X=X) #From the help-page
#I still need to consider this function (rbintrait) in more detail, to understand it better.
#Other option: y<-rbinTrait(n=1, phy=tree,beta=c(0.7),alpha=1)
table(y)
#Plot the binary states (in red/blue on the phylogeny)
plot.phylo(tree,tip.color = c("blue","red")[factor(y)])

### Including Species names:
sp <- tree$tip.label               
regre <- data.frame(sp,y)  
