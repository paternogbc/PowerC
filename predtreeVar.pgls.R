## predtreeVar.pgls: PGLS analysis accounting for variability in predictors and trees
## Author: Caterina Penone (caterina.penone@gmail.com)
## Version: 0.1
## Data created: 27.01.15
## Function under construction..!!

#Load packages
library(caper,warn.conflicts = F,quietly = T)

predtreeVar.pgls <- function(resp,pred,se.pred,tree,ntree=1,npred=1,method=c("normal","uniform"),taxa.col,lambda="ML"){
  
  # Error check
  method <- match.arg(method)
  data<-data.frame(taxa.col,resp,pred,se.pred)
  
  #Function to pick a random value in the interval
  if(method=="normal") funr <- function(a, b) {rnorm(1,a,b)}
  else
    funr <- function(a, b) {runif(1,a-b,a+b)}
  
  #Create the results data.frame
  resultados<-data.frame("n.tree"=numeric(),"n.pred"=numeric(),"estim"=numeric(),
                         "rsq"=numeric(),"pval"=numeric(),"Aicc"=numeric())
  
  #Model calculation
  counter=1
  c.data<-list()
  for (i in 1:length(ntree)) {
    for (j in 1:npred){
      
      #choose a random value in [mean-sd,mean+sd]
      data$variab<-apply(data[,c("pred","se.pred")],1,function(x)funr(x["pred"],x["se.pred"])) 
      
      #comparative data creation if tree is class=multiphylo
      if (inherits(tree, "multiPhylo")) {
        c.data[[i]]<-comparative.data(data=data, phy=tree[[i]], names.col="taxa.col", vcv=T, vcv.dim=3) ###
        ModeloSimple<- pgls(resp~variab, c.data[[i]], lambda='ML')
      }
      
      else {
        c.data<-comparative.data(data=data, phy=tree, names.col="taxa.col", vcv=T, vcv.dim=3) ###
        ModeloSimple<- pgls(resp~variab, c.data, lambda='ML')
      }
      
      #extract model coefficients
      estim<-summary(ModeloSimple)$coef[2,1]
      resq<-summary(ModeloSimple)$r.squared
      pval<-summary(ModeloSimple)$coef[2,4]
      Aicc<-ModeloSimple$aicc
      
      #write in a table
      resultados[counter,1]<- i
      resultados[counter,2]<- j    
      resultados[counter,3]<- estim
      resultados[counter,4]<- resq
      resultados[counter,5]<- pval
      resultados[counter,6]<- Aicc
      counter=counter+1
    }
  }
  
  
  #calculate mean and sd for each parameter
  #variation due to tree choice
  mean_by_tree<-aggregate(.~n.tree, data=resultados, mean)
  #variation due to continuous trait
  mean_by_randomval<-aggregate(.~n.pred, data=resultados, mean)
  
  statresults<-data.frame(mean=apply(resultados,2,mean),
                          sd_all=apply(resultados,2,sd),
                          sd_tree=apply(mean_by_tree,2,sd),
                          sd_pred=apply(mean_by_randomval,2,sd))[-(1:2),]
  
  
  output <- list(ntree=ntree,formula=formula,succesful_iterations=nrow(resultados),
                 model_results=resultados,stats=statresults)
  
  return(output)  
}