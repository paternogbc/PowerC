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

#Simulate binary trait with phylogenetic signal. #Also consider using the (primates data from corHMM) 
y<-rbinTrait(n=1, phy=tree,beta=c(0.5),alpha=1)
table(y)
#Plot the binary states (in red/blue on the phylogeny)
plot.phylo(tree,tip.color = c("blue","red")[factor(y)])

### Including Species names:
sp <- tree$tip.label               
regre <- data.frame(sp,y)  

####Some basic model running
#Run the simplest model, with 1 rate class (i.e. reduces to mk2)
system.time(mod0_hrm1 <- corHMM(phy = tree,data = regre,rate.cat = 1,node.states = "marginal",nstarts=10))
#For most uses this is realistically going to require a cluster. Here, only 50 tips and 10 restarts (more appropriate usual ~100) and takes 20 seconds for simplest model
#Other rate classes
system.time(mod0_hrm2 <- corHMM(phy = tree,data = regre,rate.cat = 2,node.states = "marginal",nstarts=10))
#Conceptually three and four are the same idea, but take very long. Consider later. 

mod0_hrm1$AICc
mod0_hrm1$solution
as.vector(mod0_hrm1$solution)
mod0_hrm1$solution.se
as.vector(mod0_hrm2$solution.se)

mod0_hrm2$AICc
mod0_hrm2$solution
as.vector(mod0_hrm2$solution)
mod0_hrm2$solution.se
as.vector(mod0_hrm2$solution.se)


########The sampling function (based on Gustavo's sampling.pgls + on my own previously written code to similar things in corHMM)

#What do we wish the function to do? 
#Within a rate class: give mean/median/se/range etc. for the AICc to check if that is somewhat stable.
#Give an estimate of the stability of the trainsition rates, and ideally something like a histogram.
#Implement an option so that people can save not only the key variables of the model, but also the actual models.
    #I think this will be interesting to many people as they will often want to actually print the triat reocnstructions and need the actual models. Particular because computatoinally intense, may be wise to store. 
#Problem: There is som many of them for high rate classes.
#Ideally: compare different rate classes
#Check histo's I made orginally for mk2 models
#Check potential for printing phylo's?

#Try to write a sampling.corHMM, comparable to what Gustavo did. 
sampling.corhmm <- function(phylogeny,data,times=20,breaks=seq(.1,.7,.1),rate.cat=1,node.states="marginal",nstarts=10)
{
  require(corHMM)
  require(ape)
  ### Basic error checking:
  if(class(phylogeny)!="phylo") stop("The phylogeny must be of the class 'phylo'. See ?corHMM")
  if(class(data)!="data.frame") stop("Data data must be of class 'data frame'. See ?corHMM")
  if(length(breaks)<2) stop("please include more than one break (eg. breaks=c(.3,.5)") 
  else
    
    # FULL MODEL calculations:
  N <- nrow(data)             # Sample size
  mod.0 <- corHMM(phy = phylogeny,data = data,rate.cat = rate.cat,node.states = node.states,nstarts=nstarts)
  AICc.0 <-mod.0$AICc
  Solution.0 <-mod.0$solution
  Solution.se.0 <-mod.0$solution
  Solution.vector.0 <-as.vector(mod.0$solution)
  solution.se.vector.0 <-as.vector(mod.0$solution)
  
  # Sampling effort analysis: 
  AICc <- as.numeric()
  Solution <- NULL
  Solution.se <- NULL
  Solution.vector <-NULL
  Solution.se.vector <-NULL
  
  # Loop:
  limit <- sort(round((breaks)*nrow(data),digits=0))
  for (i in limit){
    for (j in 1:times){ 
      exclude <- sample(1:N,i)
      crop.data <- data[-exclude,]
      small.phylogeny<-drop.tip(phy=phylogeny,tip=phylogeny$tip.label[!phylogeny$tip.label %in% crop.data$sp]) #prune tree. This will now only work if the species column is called 'sp'--> make general
      mod=try(corHMM(phy = small.phylogeny,data = crop.data,rate.cat = rate.cat,node.states = node.states,nstarts=nstarts),TRUE)
      if(isTRUE(class(mod)=="try-error")) { next } 
      else { 
        AICc <-mod$AICc
        Solution <-mod$solution
        Solution.se <-mod$solution
        Solution.vector <-as.vector(mod$solution)
        solution.se.vector <-as.vector(mod$solution)
        
        ### Storing values for each simulation
        AICc <- c(AICc,AICc)
        Solution.vector <- rbind(Solution.vector,Solution.vector)
        Solution.se.vector <- rbind(Solution.se.vector,Solution.se.vector)
              } 
            
    }
  }
    return(list(AICc.0,Solution.0,Solution.se.0,Solution.vector.0,solution.se.vector.0,
                AICc,data.frame(Solution.vector),data.frame(Solution.se.vector)))
  #Much to do here still: Implement some actual analyses of these.
  #Implement the option to save the actual models
}

#Try and run it
corhmm_samp_try<-sampling.corhmm(phylogeny = tree,data = regre,times = 10,breaks = seq(.2,.8,.2),rate.cat = 1)
corhmm_samp_try2<-sampling.corhmm(phylogeny = tree,data = regre,times = 10,breaks = seq(.2,.8,.2),rate.cat = 2)
