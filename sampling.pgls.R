## sampling.pgls: Analysis of sampling effort for PGLS linear regression
## Author: Gustavo Paterno (paternogbc@gmail.com)
## Version: 0.5
## Data created: 10.10.14
## Last update: 21.11.14

## This code is not totally checked, please be aware!

## Load required packages:
require(caper)
require(ggplot2)
require(gridExtra)

sampling.pgls <- function(formula,data,phy,times=5,breaks=c(.1,.3,.5,.7,.9),lambda=1,names.col)
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
                                        n.percent <- c(n.remov/N)
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
          result <- data.frame(intercepts,betas,p.values,n.removs,n.percent)
          ## Graphs: Estimate Betas Density plot:
          .e <- environment()
          p1 <- ggplot(result,aes(x=betas),environment=.e)+
                    geom_density(size=.5, fill="tomato",alpha=.5)+
                    xlab("Betas")+
                    geom_vline(xintercept=beta.0.low,linetype=2,color="red")+
                    geom_vline(xintercept=beta.0.up,linetype=2,color="red")+
                    geom_vline(xintercept = beta.0,color="red",linetype=2,size=2)                    
          
          ## Graphs: Estimated betas ~ % species removed
          p2 <- ggplot(result,aes(y=betas,x=n.removs))+
                    geom_point(size=3,alpha=.7)+
                    ylab("Estimated Beta")+
                    xlab("% of Species Removed ")+
                    geom_hline(yintercept=beta.0.low,linetype=2,color="red")+
                    geom_hline(yintercept=beta.0.up,linetype=2,color="red")+
                    geom_hline(yintercept=beta.0,linetype=2,color="red",size=1.1)
          
          ## Mean estimated Betas:
          med <- with(result,tapply(betas,n.removs,mean))
          Sdev <- with(result,tapply(betas,n.removs,sd))
          n.sp <- nrow(c.data$data)-as.numeric(rownames(med))
          result.med <- data.frame(med,Sdev,n.sp,beta.0)
          p3 <- ggplot(result.med,aes(y=med,x=n.sp),environment=.e)+
                    geom_point(size=3,alpha=.7)+
                    scale_x_continuous(breaks=c(n.sp,nrow(c.data$data)))+
                    geom_errorbar(aes(ymin=med-Sdev, ymax=med+Sdev), width=.1)+
                    geom_hline(yintercept= beta.0,linetype=2,color="red")+
                    xlab("Number of Species") + ylab("Estimated Beta (+-SD)")+
                    geom_point(aes(x=nrow(c.data$data),y=beta.0,size=3,colour="red"))+
                    theme(legend.position="none")
                      
          ## Power Analysis:
          simu.sig <- result$p.values > .05
          result$simu.sig <- simu.sig
          power <- 1-(with(result,tapply(simu.sig,n.removs,sum)))/times
          power.tab <- data.frame(breaks,power)
          p4 <-ggplot(power.tab,aes(y=power,x=breaks))+
                    scale_y_continuous(limits=c(0,1))+ 
                    scale_x_continuous(breaks=breaks)+          
                    xlab("% Species removed")+
                    geom_point(size=5,colour="red")+
                    geom_line(colour="red")
          grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
          return(result)
}

