PowerC
======

PowerC is a package of R functions to perform power analysis in comparative methods

###Functions:

  power.pgls - Simple power analysis of PGLS linear regression

###Description:

	This function performs a simple power analysis of a linear regression in PGLS. The function compare the full model (with all data samples) with subsets of the full dataset. In each run, a number of species are removed at random. Species are removed in intervals ranging from 5%-90% of the number of species in the comparatve dataset). The user can specify the number of replicates per interval. Estimates and significance of each simulation is stored for later comparison.

###Usage:
	
	power.pgls(x,y,data,phy,times=5,breaks=c(.1,.3,.5,.7,.9),lambda="ML",names.col)

###Arguments:

	formula		A model formula
	data 		A data frame containing data for each tip on the phylogeny
	phy         A phylogeny (class 'phylo') to be matched to the data above.
	names.col   The name of a column in the provided data frame that will be used to match data rows to phylogeny tips.
	lambda		A value for the lambda transformation.
	breaks		Numeric vector with percentages of data to be removed. (eg. breaks= c(.1,.5,.9) will run simulations with 10%,50% and 90% of species removed	
	times  		Number of simulations for each break
		

###Details:
	
	This function only works for simple linear regression (Y ~ bx + a). Future implementations will deal with more complex models.

###Value:

	The ‘power.pgls’ function returns an object of class data.frame containing the following:


	result : 	A data.frame with the estimates (intercept and beta), p-values, and number of species removed in each simulation


 
Warning: This code still needs to be fully checked, be aware.


Note:


###Author(s):

	Gustavo Brant Paterno
	paternogbc@gmail.com

###References:

	R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic analysis and comparative data: A test and review of evidence. American Naturalist, 160:712-726, 2002.


###See Also:

	pgls, comparative.data


### Exampless
	### Required packages:
	library(caper);library(phytools);library("phylolm")
	set.seed(111)
	### Simulating tree
	tree<-pbtree(n=100)    
	### Simulating response variable with phylogenetic signal
	y <- rTrait(n=1, tree, model=c("lambda"),parameters=list(lambda=.8))  
	### Simulating explanatory variable
	x <- varY + rnorm(100,mean(varY),1)     
	### Including Species names
	sp <- tree$tip.label               
	data <- data.frame(sp,y,x,x2)   

	### Organizing comparative data for pgls:
	c.data <- comparative.data(data=data,phy=tree,vcv=T,vcv.dim=3,names.col="sp")

	### Fiting lm and pgls regressions:
	mod.pgls <- pgls(y ~ x, c.data,"ML")
	summary(mod.pgls)
  
	### Regression Plots:
	plot(y ~ x,pch=16)
	abline(mod.pgls)
     
	### Perform sensitive analysis:
	result <- power.pgls(y~x,data=data,phy=tree,times=10,breaks=seq(.1,.9,.1),lambda="ML",names.col="sp")
	result2 <- power.pgls(y~x,data=data,phy=tree,times=50,breaks=c(.1,.5,.9),lambda="ML",names.col="sp")

	### Visualizing results:
	med <- with(result,tapply(betas,n.removs,mean))
	Sdev <- with(result,tapply(betas,n.removs,sd))

	n.re <- as.numeric(rownames(med))
	n.sp <- nrow(c.data$data)-n.re 
	result.med <- data.frame(med,Sdev,n.sp,n.re)
	beta.mod.pgls <- summary(mod.pgls)[[c(5,2)]] 

	ggplot(result.med,aes(y=med,x=n.sp))+
          geom_point(size=3,alpha=.7)+
          scale_x_continuous(breaks=seq(5, nrow(c.data$data), 5))+
          geom_errorbar(aes(ymin=med-Sdev, ymax=med+Sdev), width=.1)+
          geom_hline(yintercept= beta.mod.pgls,linetype=2,color="red")+
          xlab("Number of Species") + ylab("Estimated Beta")+
          geom_point(aes(x=nrow(c.data$data),y=beta.mod.pgls,size=3,colour="red"))+
          theme(legend.position="none",
                axis.text=element_text(size=14),       
                axis.title=element_text(size=20))


	

