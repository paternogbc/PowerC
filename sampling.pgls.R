## sampling.pgls: Analysis of sampling effort for PGLS linear regression
## Author: Gustavo Paterno (paternogbc@gmail.com)
## Version: 0.7
## Data created: 10.12.14

## This code is not totally checked, please be aware!

## Load required packages:
library(caper)

sampling.pgls <- function(formula,data,phy,times=20,breaks=c(.1,.3,.5,.7,.9),lambda="ML",names.col)
{
          ### Basic error checking:
          if(class(formula)!="formula") stop("Please formula must be class 'forumla'")
          if(class(phy)!="phylo") stop("phy object must be of class 'phylo'")
          if(class(data)!="data.frame") stop("data data must be of class 'data.frame'")
          if(length(breaks)<2) stop("please include more then one break (eg. breaks=c(.3,.5)") 
          else
                    
          # FULL MODEL calculations:
          c.data <- comparative.data(phy=phy,data=data,names.col=sp,vcv=T,vcv.dim=3)
          N <- nrow(c.data$data)             # Sample size
          mod.0 <- pgls(formula, data=c.data,lambda=lambda)
          sumMod <- summary(mod.0)
          intercept.0 <-    sumMod[[c(5,1)]] # Intercept (full model)
          beta.0 <-    sumMod[[c(5,2)]]      # Beta (full model)
          pval.0 <-    sumMod[[c(5,8)]]      # p.value (full model)
          sd.beta.0 <- sumMod[[c(5,4)]]      # Standart Error (full model)
          df.0 <- sumMod[[2]][2] # Degree if Freedon (full model))
          beta.IC <- qt(0.975,df.0)*sd.beta.0 # Beta CI (full model)
          beta.0.low <- beta.0 - beta.IC  # Low limit of beta CI (full model)
          beta.0.up <- beta.0 + beta.IC   # Up limit of beta CI (full model)
          
          # Sampling effort analysis: 
          intercepts <- as.numeric()
          betas <- as.numeric()
          p.values <- as.numeric()
          n.removs <- as.numeric()
          n.percents <- as.numeric()
          
          # Loop:
          limit <- sort(round((breaks)*nrow(c.data$data),digits=0))
          for (i in limit){
                    for (j in 1:times){ 
                              out <- sample(1:N,i)
                              crop.data <- comparative.data(phy=phy,data=data[-out,],names.col=sp,vcv=T,vcv.dim=3)
                              mod=try(pgls(formula, data=crop.data,lambda),TRUE)
                              if(isTRUE(class(mod)=="try-error")) { next } 
                              else { 
                                        ### Calculating model estimates:
                                        sum.Mod <- summary(mod)
                                        beta <-    sum.Mod[[c(5,2)]]     # Beta
                                        intercept <-    sum.Mod[[c(5,1)]]# Intercept
                                        pval <-    sum.Mod[[c(5,8)]] # p.value
                                        n.remov <- i
                                        n.percent <- n.remov/N
                                        rep <- j
                                        
                                        ### Storing values for each simulation
                                        intercepts <- c(intercepts,intercept)
                                        betas <- c(betas,beta)
                                        p.values <- c( p.values,pval)
                                        n.removs <- c(n.removs,n.remov)
                                        n.percents <- c(n.percents,n.percent)
                              } 
                              
                              
                    }
          }
          
          
          # Data frame with results:
          estimates <- data.frame(intercepts,betas,p.values,n.removs,n.percents)
          
          ## Power Analysis:
          times <- table(estimates$n.removs)[1]
          breaks <- unique(estimates$n.percents)
          simu.sig <- estimates$p.values > .05
          estimates$simu.sig <- simu.sig
          power <- 1-(with(estimates,tapply(simu.sig,n.removs,sum)))/times
          power.tab <- data.frame(percent_sp_removed=breaks,power)
          estimates <- estimates[,-ncol(estimates)]
          
          param0 <- data.frame(beta.0,intercept.0)
          beta_IC <- data.frame(beta.low=beta.0.low,beta.up=beta.0.up)
          return(list(original_model_estimates=param0,original_beta_95_IC=beta_IC,results=estimates,power_analysis=power.tab))

  
}

