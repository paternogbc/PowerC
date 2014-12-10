## influential.pgls: Influential analysis for PGLS linear regression
## Author: Gustavo Paterno (paternogbc@gmail.com)
## Version: 0.1
## Data created: 26.11.14
## Last update: 

## This code is not totally checked, please be aware!

## Load required packages:
require(caper)
require(ggplot2)
require(gridExtra)

influence.pgls <- function(formula,data,phy,lambda=1,names.col)
{
          ### Basic error checking:
          if(class(formula)!="formula") stop("Please formula must be class 'forumla'")
          if(class(phy)!="phylo") stop("phy object must be of class 'phylo'")
          if(class(data)!="data.frame") stop("data data must be of class 'data.frame'")
          else
                    
          # FULL MODEL calculations:
          c.data <- comparative.data(phy=phy,data=data,names.col=sp,vcv=T,vcv.dim=3)
          N <- nrow(c.data$data)             # Sample size
          mod.0 <- pgls(formula, data=c.data,lambda=lambda)
          sumMod <- summary(mod.0)
          intercept.0 <-    sumMod[[c(5,1)]] # Intercept (full model)
          beta.0 <-    sumMod[[c(5,2)]]      # Beta (full model)
          pval.0 <-    sumMod[[c(5,8)]]      # p.value (full model)
          
          # Sampling effort analysis:
          betas <- as.numeric()
          intercepts <- as.numeric()
          DFbetas <- as.numeric()
          DFintercepts <- as.numeric()
          p.values <- as.numeric()
          species <- as.character()
          
          # Loop:
          for (i in 1:nrow(c.data$data)){
                    crop.data <- comparative.data(phy=phy,data=data[-i,],names.col=sp,vcv=T,vcv.dim=3)
                    mod=try(pgls(formula, data=crop.data,lambda),TRUE)
                    if(isTRUE(class(mod)=="try-error")) { i = i-1 
                                                          next } 
                    else { 
                                        ### Calculating model estimates:
                                        sum.Mod <- summary(mod)
                                        beta <-    sum.Mod[[c(5,2)]]     # Beta
                                        intercept <-    sum.Mod[[c(5,1)]]# Intercept
                                        pval <-    sum.Mod[[c(5,8)]] # p.value
                                        DFbeta <- beta - beta.0
                                        DFint  <- intercept - intercept.0
                                        sp <- phy$tip.label[i]
                                        ### Storing values for each simulation
                                        betas <- c(betas,beta)
                                        intercepts <- c(intercepts,intercept)
                                        DFbetas <- c(DFbetas,DFbeta)
                                        DFintercepts <- c(DFintercepts,DFint)
                                        species <- c(species,sp)
                                        p.values <- c( p.values,pval)
                                        
                              } 
          }
          # Data frame with results:
          estimates <- data.frame(species,betas,DFbetas,intercepts,DFintercepts,p.values) 
          param <- data.frame(intercept.0,beta.0)
          influ.sp <- as.character(estimates[order(estimates$DFbetas,decreasing=T)[1:5],]$species)
          
          return(list(estimates=estimates,ori_coeff=param,data=data,formula=formula,influ.sp=influ.sp))
          
}

