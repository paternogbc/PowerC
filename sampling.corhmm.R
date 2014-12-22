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
#Implement: 
    #Calculate mean/se/median of trans rates per break --> also (or only?) in visualising funciton
    #Idem for AIC
    #dplyr will be useful here. 
#Ideally: compare different rate classes
  #For instance, iterate over different rate classes and see what happens?
  #Try if this is possible (it should be, but time consuming)

#For visualisation
  #Check histo's I made orginally for mk2 models
  #Check potential for printing phylo's? --> Write a seperate function to iterate over stored models and print the phylogenies. 

#Implement an option for mk2 in diversitree as well. 
  #Check other main packages for binary trait reconstruction (i.e. ape etc.).

  
#Try to write a sampling.corHMM, comparable to what Gustavo did. 
sampling.corhmm <- function(phylogeny,data,times=10,breaks=seq(.1,.9,.1),rate.cat=1,node.states="marginal",nstarts=10,save.models=F)
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
  Solution.se.0 <-mod.0$solution.se
  Solution.vector.0 <-as.vector(mod.0$solution)
  Solution.se.vector.0 <-as.vector(mod.0$solution.se)
  
  # Sampling effort analysis: 
  AICc <- as.numeric()
  Solution.vector <-NULL
  Solution.se.vector <-NULL
  n.removs <- as.numeric()
  n.percents <- as.numeric()
  mods<-list()
  
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
        AICc_new <-mod$AICc
        Solution.vector_new <-as.vector(mod$solution)
        Solution.se.vector_new <-as.vector(mod$solution.se)
        n.remov <- i
        n.percent <- n.remov/N
        
        ### Storing values for each simulation
        AICc <- c(AICc,AICc_new)
        Solution.vector <- rbind(Solution.vector,Solution.vector_new)
        Solution.se.vector <- rbind(Solution.se.vector,Solution.se.vector_new)
        n.removs <- c(n.removs,n.remov)
        n.percents <- c(n.percents,n.percent)
        if(isTRUE(save.models)) {mods<-c(mods,list(mod))} 
              } 
            
    }
  }
    #Data frame with results
  estimates<-data.frame(AICc,Solution.vector,Solution.se.vector,n.removs,n.percents)
  original_trans_vector<-data.frame(t(Solution.vector.0),t(Solution.se.vector.0))
  
 return(list(AICc_original=AICc.0,Transition_rates_original=Solution.0,
              SE_transition_rates_original=Solution.se.0,
              Original_rates_vectorised=original_trans_vector,
              results=estimates,models=mods))
  
 }

#Try and run it (#Low replicate numbers, only one rate class for simplicity)
corhmm_samp_try<-sampling.corhmm(phylogeny = tree,data = regre,times = 2,breaks = seq(.2,.8,.2),rate.cat = 1,save.models = T)
#With more rate classes
#corhmm_samp_try2<-sampling.corhmm(phylogeny = tree,data = regre,times = 3,breaks = seq(.2,.8,.2),rate.cat = 2)
