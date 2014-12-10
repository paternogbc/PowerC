## plot.power.pgls: Plot influential and sampling effort analysis for pgls
## Author: Gustavo Paterno (paternogbc@gmail.com)
## Version: 0.1
## Data created: 10.10.14
## Last update: 21.11.14

## This code is not totally checked, please be aware!

require(ggplot2)
require(gridExtra)

### Start:
plot.power.pgls <- function(x, method="sampling"){
          if (method == "sampling"){
                    result <- x[[1]]
                    beta.0 <- as.numeric(x[[2]][1])
                    beta.0.low <- as.numeric(x[[2]][2])
                    beta.0.up <- as.numeric(x[[2]][3])
                    .e <- environment()
                    p1 <- ggplot(result,aes(x=betas),environment=.e)+
                              geom_density(size=.5, fill="tomato",alpha=.5)+
                              xlab("Betas")+
                              geom_vline(xintercept=beta.0.low,linetype=2,color="red")+
                              geom_vline(xintercept=beta.0.up,linetype=2,color="red")+
                              geom_vline(xintercept = beta.0,color="red",linetype=2,size=2)                    
                    
                    ## Graphs: Estimated betas ~ % species removed
                    p2 <- ggplot(result,aes(y=betas,x=n.percents))+
                              geom_point(size=3,alpha=.7)+
                              scale_x_continuous(breaks=result$n.percents)+          
                              ylab("Estimated Betas")+
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
                              xlab("Number of Species") + ylab("Mean Estimated Betas (+-SD)")+
                              geom_point(aes(x=nrow(c.data$data),y=beta.0,size=3,colour="red"))+
                              theme(legend.position="none")
                    
                    ## Power Analysis:
                    times <- table(result$n.removs)[1]
                    breaks <- unique(result$n.percents)
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
          }
          if (method == "influential"){
                    result <- x[[1]]
                    intercept.0 <-  as.numeric(x[[2]][1])
                    beta.0 <-  as.numeric(x[[2]][2])
                    p1 <- ggplot(result,aes(x=betas),environment=.e)+
                              geom_density(size=.5, fill="tomato",alpha=.5)+
                              xlab("Betas")+
                              geom_vline(xintercept = beta.0,color="red",linetype=2,size=2)
                    p2 <- ggplot(result,aes(x=intercepts),environment=.e)+
                              geom_density(size=.5, fill="tomato",alpha=.5)+
                              xlab("Intercepts")+
                              geom_vline(xintercept = intercept.0,color="red",linetype=2,size=2)
                    p3<-qplot(abs(result$betas - beta.0),xlab=("DF Betas"))
                    p4<-qplot(abs(result$intercepts - intercept.0),xlab=("DF Intercepts"))
                    result.tab <- data.frame(result,x$data[all.vars(x$formula)])
                    
                    p5<-ggplot(result.tab,aes(y=y,x=x,colour=abs(DFbetas)))+
                              geom_point(size=3)+
                              scale_colour_gradient( low="black", high="red", space="Lab",name="DF Betas")
                    p6<-ggplot(result.tab,aes(y=y,x=x,colour=abs(DFintercepts)))+
                              geom_point(size=3)+
                              scale_colour_gradient( low="black", high="red", space="Lab",name="DF Intercepts")          
                              grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
          }
                    
                    
          
}
          
