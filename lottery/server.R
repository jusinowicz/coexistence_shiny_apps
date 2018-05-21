library(shiny)
########################
#
# This is R code to run :
# 	 Lottery model with stochastic growth rates
# 
# This code always assumes 2 species and that space is proportionate
# (i.e. is constrained between 0 and 1) 
#
# 1. The lottery model was first used by Chesson and Warner (1981) 
#    to model recruitment in coral-reef fish. Competition is for
#    open space (made available through adult mortality). It is
#    based on the proportion of potential recruits that each species
#    has produced in a generation. 

# 2. Stochastic reproduction can be modelled in a relatively 
#    straightforward way.

#Load library MASS For multivariate random variables
library(MASS)
###################################
# 1. Leslie Gower model with fluctuating reproduction
######################################################

#This produces a column vector for each species. 
#That is, the dimensions are such that rows = time. 
pop_lott = function (ia,ngens, rs, surv,cor,var) {
	
	#par(mar = c(0,0,0,0))

	ia=as.matrix(c(ia,1-ia) )
	rs=as.matrix(rs)
	surv=as.matrix(surv)
	var=as.matrix(var)


	#Initialize populations
	lott1 = matrix(0, ngens, 1)
	lott2 = lott1

	#Demographic parameters. 
	#Reproduction rates are multivariate random variables now! 
	#Means
	m1=rs[1]
	m2=rs[2]
	means12=c(m1,m2)

	#Variances
	v1=var[1]
	v2=var[2]
	vs12=c(v1,v2)

	#correlation/covariance
	rho12=cor
	cov12 = matrix(c(v1,(rho12*sqrt(v1)*sqrt(v2)),(rho12*sqrt(v2)*sqrt(v1)),v2),2,2) 

	#Transform specified means, variance to LOGNORMAL means, variances
	MLog=log((means12^2)/sqrt(vs12+means12^2))
	VLog=log(vs12/(means12^2)+1);
	covLog=log(1+(cov12)/abs(means12%*%t(means12)));

	#Now use mvrnorm to generate the entire time series of reproduction. 
	#repro12 = mvrnorm (ngens, means12, cov12)
	repro12 = mvrnorm (ngens, MLog, covLog)
	repro12log=exp(repro12)
	r1= repro12log[,1]
	r2= repro12log[,2]

	#Survival rates
	s1 = surv[1]
	s2 = surv[2]

	#Important: 
	#Invasion is simulated in a single loop. The resident
	#and invader designation is determined by the initial conditions (ICs). 

	#Give each species its initial abundance
	lott1[1] = ia[1]
	lott2[1] = ia[2]


	#Population growth
	for ( t in 1:ngens) {

		lott1[t+1] = s1*lott1[t]+(1-sum(surv*c(lott1[t],lott2[t])))*r1[t]*lott1[t]/(r1[t]*lott1[t]+r2[t]*lott2[t])
		lott2[t+1] = s2*lott2[t]+(1-sum(surv*c(lott1[t],lott2[t])))*r2[t]*lott2[t]/(r1[t]*lott1[t]+r2[t]*lott2[t])

	}

	#Visualize,
	#jpeg( filename="lgnr1f_both_co.jpg", width=5, height=5, units = "in", pointsize=12, quality=100,res=300)
	plot(lott2, t="l",col="red",xlab="Generations", ylab="Population size",ylim=c(0,max(cbind(lott1,lott2))))
	lines(lott1, col="green")

#dev.off
}


######################################################
# 1. Leslie Gower model with fluctuating reproduction
######################################################
#
# Animate the population growth in space: 
#

pop_space_lott= function (ia,ngens, rs, surv,cor,var) {
	

	#For spatial visualization
	np=100

	#par(mar = c(0,0,0,0))
	
	ia=as.matrix(c(ia,1-ia) )
	rs=as.matrix(rs)
	surv=as.matrix(surv)
	var=as.matrix(var)


	#Initialize populations
	lott1 = matrix(0, ngens, 1)
	lott2 = lott1

	#Demographic parameters. 
	#Reproduction rates are multivariate random variables now! 
	#Means
	m1=rs[1]
	m2=rs[2]
	means12=c(m1,m2)

	#Variances
	v1=var[1]
	v2=var[2]
	vs12=c(v1,v2)

	#correlation/covariance
	rho12=cor
	cov12 = matrix(c(v1,(rho12*sqrt(v1)*sqrt(v2)),(rho12*sqrt(v2)*sqrt(v1)),v2),2,2) 

	#Transform specified means, variance to LOGNORMAL means, variances
	MLog=log((means12^2)/sqrt(vs12+means12^2))
	VLog=log(vs12/(means12^2)+1);
	covLog=log(1+(cov12)/abs(means12%*%t(means12)));

	#Now use mvrnorm to generate the entire time series of reproduction. 
	#repro12 = mvrnorm (ngens, means12, cov12)
	repro12 = mvrnorm (ngens, MLog, covLog)
	repro12log=exp(repro12)
	r1= repro12log[,1]
	r2= repro12log[,2]

	#Survival rates
	s1 = surv[1]
	s2 = surv[2]

	#Important: 
	#Invasion is simulated in a single loop. The resident
	#and invader designation is determined by the initial conditions (ICs). 

	#Give each species its initial abundance
	lott1[1] = ia[1]
	lott2[1] = ia[2]


	#Population growth
	for ( t in 1:ngens) {

		#Spread population 1 out in space: 
		prop1 = lott1[t]/(lott1[t]+lott2[t]) #Proportion of space owned by 1
		pop_space = matrix( as.numeric( runif(np^2) < prop1),np,np) #Randomly assign cells to 1 based on proportion

		#Plot and delay each new plot to create appearance of animation

		image(pop_space)
		#date_time<-Sys.time()
		#while((as.numeric(Sys.time()) - as.numeric(date_time))<0.25){} #dummy while loop

		#Populations for next time step
		lott1[t+1] = s1*lott1[t]+(1-sum(surv*c(lott1[t],lott2[t])))*r1[t]*lott1[t]/(r1[t]*lott1[t]+r2[t]*lott2[t])
		lott2[t+1] = s2*lott2[t]+(1-sum(surv*c(lott1[t],lott2[t])))*r2[t]*lott2[t]/(r1[t]*lott1[t]+r2[t]*lott2[t])


	}

	return(cbind(lott1,lott2)) 

}


#################################################
# Shiny server logic
#################################################
# This section takes inputs from the ui and displays the desired
# information using the code from Saavedra et al. above

shinyServer( function(input, output) {
	
	output$pop_lott = renderPlot({
		pop_lott(input$ia, input$ngens, input$rs, input$surv,input$cor, input$var )
		})

	output$pop_space_lott = renderPlot({ 
		pop_space_lott (input$ia, input$ngens, input$rs, input$surv,input$cor, input$var )
		})

	#Download the population data
	datasetInput = reactive( pop_space_lott (input$ia, input$ngens, input$rs, input$surv,input$cor, input$var ))

	output$downloadData =	 downloadHandler(
    	filename= "popLott_data.csv",
    	
    	content = function(file) {
      	write.csv(datasetInput(), file, row.names = FALSE)
    	})

})