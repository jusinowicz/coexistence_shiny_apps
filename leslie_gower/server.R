library(shiny)
########################
#
# This is R code to run :
# 	 Leslie-Gower competition model with stochastic growth rates
# 
# This code always assumes 2 species 
#
# 1. The Leslie-Gower model is also sometimes called the multispecies 
#    Beverton-Holt model. It has been used by Jonathan and others as 
#    a model of annual plants.
#  	lg1, lg2: populations of species 1,2 
#	rlg1, rlg2: reproduction rates	
#       alpha_ij: competition coefficients, where 
#		 the notation "ij" gives the competitive
#                effect of species j on i.  
#
# 2. Stochastic reproduction can be modelled in a relatively 
#    straightforward way.

#Load library MASS For multivariate random variables
library(MASS)
###################################
# 1. Leslie Gower model with fluctuating reproduction
######################################################

#This produces a column vector for each species. 
#That is, the dimensions are such that rows = time. 
pop_lgf = function (ia,ngens, rs, alphamat,surv,cor,var) {
	
	#par(mar = c(0,0,0,0))

	ia=as.matrix(ia)
	rs=as.matrix(rs)
	alphamat=as.matrix(alphamat)
	surv=as.matrix(surv)
	var=as.matrix(var)


	#Initialize populations
	lg1f = matrix(0, ngens, 1)
	lg2f = lg1f

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
	r1lgf= repro12log[,1]
	r2lgf= repro12log[,2]

	#Survival rates
	s1f = surv[1]
	s2f = surv[2]

	#Intra/Inter specific competition
	alpha_11 = alphamat[1,1]
	alpha_22 = alphamat[2,2]
	alpha_21 = alphamat[2,1]
	alpha_12 = alphamat[1,2]

	#Important: 
	#Invasion is simulated in a single loop. The resident
	#and invader designation is determined by the initial conditions (ICs). 

	#Give each species its initial abundance
	lg1f[1] = ia[1]
	lg2f[1] = ia[2]


	#Population growth
	for ( t in 1:ngens) {

		lg1f[t+1] = s1f*lg1f[t]+r1lgf[t]*lg1f[t]/(1+lg1f[t]*alpha_11+lg2f[t]*alpha_12)
		lg2f[t+1] = s2f*lg2f[t]+r2lgf[t]*lg2f[t]/(1+lg2f[t]*alpha_22+lg1f[t]*alpha_21)

	}

	#Visualize,
	#jpeg( filename="lgnr1f_both_co.jpg", width=5, height=5, units = "in", pointsize=12, quality=100,res=300)
	plot(lg2f, t="l",col="red",xlab="Generations", ylab="Population size",ylim=c(0,max(cbind(lg1f,lg2f))))
	lines(lg1f, col="green")

#dev.off
}


######################################################
# 1. Leslie Gower model with fluctuating reproduction
######################################################
#
# Animate the population growth in space: 
#

pop_space_lgf = function (ia,ngens, rs, alphamat,surv,cor,var) {
	
	#par(mar = c(0,0,0,0))
	
	ia=as.matrix(ia)
	rs=as.matrix(rs)
	alphamat=as.matrix(alphamat)
	surv=as.matrix(surv)
	var=as.matrix(var)

	#For spatial visualization
	np=100

	#Initialize populations
	lg1f = matrix(0, ngens, 1)
	lg2f = lg1f

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
	r1lgf= repro12log[,1]
	r2lgf= repro12log[,2]

	#Survival rates
	s1f = surv[1]
	s2f = surv[2]

	#Intra/Inter specific competition
	alpha_11 = alphamat[1,1]
	alpha_22 = alphamat[2,2]
	alpha_21 = alphamat[2,1]
	alpha_12 = alphamat[1,2]

	#Important: 
	#Invasion is simulated in a single loop. The resident
	#and invader designation is determined by the initial conditions (ICs). 

	#Give each species its initial abundance
	lg1f[1] = ia[1]
	lg2f[1] = ia[2]


	#Population growth
	for ( t in 1:ngens) {

		#Spread population 1 out in space: 
		prop1 = lg1f[t]/(lg1f[t]+lg2f[t]) #Proportion of space owned by 1
		pop_space = matrix( as.numeric( runif(np^2) < prop1),np,np) #Randomly assign cells to 1 based on proportion

		#Plot and delay each new plot to create appearance of animation

		image(pop_space)
		#date_time<-Sys.time()
		#while((as.numeric(Sys.time()) - as.numeric(date_time))<0.25){} #dummy while loop

		#Populations for next time step
		lg1f[t+1] = s1f*lg1f[t]+r1lgf[t]*lg1f[t]*1/(1+lg1f[t]*alpha_11+lg2f[t]*alpha_12)
		lg2f[t+1] = s2f*lg2f[t]+r2lgf[t]*lg2f[t]*1/(1+lg2f[t]*alpha_22+lg1f[t]*alpha_21)


	}

	return(cbind(lg1f,lg2f)) 

}


#################################################
# Shiny server logic
#################################################
# This section takes inputs from the ui and displays the desired
# information using the code from Saavedra et al. above

shinyServer( function(input, output) {
	
	output$pop_lgf = renderPlot({
		pop_lgf(input$ia, input$ngens, input$rs, input$alphamat,input$surv,input$cor, input$var )
		})

	output$pop_space_lgf = renderPlot({ 
		pop_space_lgf (input$ia, input$ngens, input$rs, input$alphamat,input$surv,input$cor, input$var )
		})

	#Download the population data
	datasetInput = reactive( pop_space_lgf (input$ia, input$ngens, input$rs, input$alphamat,input$surv,input$cor, input$var ))

	output$downloadData =	 downloadHandler(
    	filename= "popLG_data.csv",
    	
    	content = function(file) {
      	write.csv(datasetInput(), file, row.names = FALSE)
    	})

})