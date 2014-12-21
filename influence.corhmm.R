#Gijsbert Werner, Vrije Universiteit Amsterdam
#Email: g.d.a.werner@vu.nl

#December 2014

#Create a R-function to analyse sensitivity to data uncertainty in binary trait evolution (specifically in corHMM)
library(phytools);library(phylolm);library(corHMM)

#Simulate a tree
set.seed(111)
N <- 50 # Number of species
### Simulating tree
tree<-pbtree(n=N)    
plot(tree)

#Everything above this, doesn't need to be redone, it's already in Gustavo's code
#Simulate binary trait with phylogenetic signal. #Also consider using the (primates data from corHMM) 
y<-rbinTrait(n=1, phy=tree,beta=c(0.5),alpha=1)
table(y)
#Plot the binary states (in red/blue on the phylogeny)
plot.phylo(tree,tip.color = c("blue","red")[factor(y)])

### Including Species names:
sp <- tree$tip.label               
regre <- data.frame(sp,y)  

#Run the simplest model, with 1 rate class (i.e. reduces to mk2)
system.time(mod0_hrm1 <- corHMM(phy = tree,data = regre,rate.cat = 1,node.states = "marginal",nstarts=10))
#For most uses this is realistically going to require a cluster. Here, only 50 tips and 10 restarts (more appropriate usual ~100) and takes 20 seconds for simplest model
#Other rate classes
system.time(mod0_hrm2 <- corHMM(phy = tree,data = regre,rate.cat = 2,node.states = "marginal",nstarts=10))
system.time(mod0_hrm3 <- corHMM(phy = tree,data = regre,rate.cat = 3,node.states = "marginal",nstarts=10))
system.time(mod0_hrm4 <- corHMM(phy = tree,data = regre,rate.cat = 4,node.states = "marginal",nstarts=10))

