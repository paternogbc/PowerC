## power.pgls: Simple power analysis of PGLS linear regression
## Author: Gustavo Paterno (paternogbc@gmail.com)
## Version: 0.5
## Data created: 10.10.14
## Last update: 21.11.14

## This code is not totally checked, please be careful!

## Load required packages:
require(caper)
require(ggplot2)
require(gridExtra)

power.pgls <- function(formula,data,phy,times=5,breaks=c(.1,.3,.5,.7,.9),lambda=1,names.col)
{
          ### Error checking:
          if(class(phy)!="phylo") stop("phy object must be of class 'phylo'")
          if(class(data)!="data.frame") stop("data data must be of class 'data.frame'")
          if(length(breaks)<2) stop("please include more then one break (eg. breaks=c(.3,.5)") 
          else
                    
                    formula <- as.formula(formula)
          # FULL MODEL calculations:
          c.data <- comparative.data(phy=phy,data=data,names.col=sp,vcv=T,vcv.dim=3)
          mod.0 <- pgls(formula, data=c.data,lambda=lambda)
          intercept.0 <-    summary(mod.0)[[c(5,1)]] # Intercept (full model)
          beta.0 <-    summary(mod.0)[[c(5,2)]]      # Beta (full model)
          pval.0 <-    summary(mod.0)[[c(5,8)]] # p.value (full model)
          sd.beta.0 <- summary(mod.0)[[c(5,4)]] # Standart Error (full model)
          lambda.0 <-      mod.0$param[2]       # Lambda estimate (full model)
          lambda.0.CI <- summary(mod.0)$param.CI[2]$lambda$ci.val # Lambda CI
          lambda.0.CI[is.na(lambda.0.CI)] <- 0
          df.0 <- summary(mod.0) [[2]][2] # Degree if Freedon (full model))
          N <- nrow(c.data$data)          # Sample size
          beta.IC <- qt(0.975,df.0)*sd.beta.0 # Beta CI (full model)
          beta.0.low <- beta.0 - beta.IC  # Low limit of beta CI (full model)
          beta.0.up <- beta.0 + beta.IC   # Up limit of beta CI (full model)
          
          # Sensitive Analysis 
          # Variables to be stored:
          intercepts <- as.numeric()
          betas <- as.numeric()
          p.values <- as.numeric()
          n.removs <- as.numeric()
          replicates <- as.numeric()
          
          # Loops:
          limit <- sort(round((1-breaks)*nrow(c.data$data),digits=0))
          
          for (i in limit){
                    for (j in 1:times){ 
                              out <- sample(1:N,i)
                              crop.data <- comparative.data(phy=phy,data=data[-out,],names.col=sp,vcv=T,vcv.dim=3)
                              mod=try(pgls(formula, data=crop.data,lambda),TRUE)
                              if(isTRUE(class(mod)=="try-error")) { next } 
                              else { 
                                        ### Calculating model estimates:
                                        beta <-    summary(mod)[[c(5,2)]]     # Beta
                                        intercept <-    summary(mod)[[c(5,1)]]# Intercept
                                        pval <-    summary(mod)[[c(5,8)]] # p.value
                                        n.remov <- i
                                        rep <- j
                                        
                                        ### Storing values for each simulation
                                        intercepts <- c(intercepts,intercept)
                                        betas <- c(betas,beta)
                                        p.values <- c( p.values,pval)
                                        n.removs <- c(n.removs,n.remov)
                                        replicates <- c(replicates,rep)
                              } 
                              
                              
                    }
          }
          
          
          # Data frame with results:
          result <- data.frame(intercepts,betas,p.values,n.removs,replicates)
          ## Graphs: Estimate Betas Density plot:
          .e <- environment()
          m <- ggplot(result,aes(x=betas),environment=.e)+
                    geom_density(size=.5, fill="tomato",alpha=.5)+
                    xlab("Betas")+
                    geom_vline(xintercept=beta.0.low,linetype=2,color="red")+
                    geom_vline(xintercept=beta.0.up,linetype=2,color="red")+
                    geom_vline(xintercept = beta.0,color="red",linetype=2,size=2)                    
          
          ## Graphs: Estimated betas ~ % species removed
          g <- ggplot(result,aes(y=betas,x=n.removs))+
                    geom_point(size=3,alpha=.7)+
                    ylab("Estimated Beta")+
                    xlab("% of Species Removed ")+
                    geom_hline(yintercept=beta.0.low,linetype=2,color="red")+
                    geom_hline(yintercept=beta.0.up,linetype=2,color="red")+
                    geom_hline(yintercept=beta.0,linetype=2,color="red",size=1.1)
          grid.arrange(m,g,ncol=2)
}

