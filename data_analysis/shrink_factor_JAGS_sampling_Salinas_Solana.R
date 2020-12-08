
if(!require("runjags")){
  install.packages("runjags")
  library(runjags)}

internodes <- read.csv("internode_length_rhizome_and_lateral.csv",header=T)

data_list <- list(
  individual=internodes$Individual,
  ninds =length(unique(internodes$Individual)),
  type = internodes$rhizome_0_else_lateral_1,
  length = internodes$internode_length,
  ntotal = dim(internodes)[1]) 

model_shrink <-"model{
for (i in 1:ntotal){
length[i]~dnorm(mu[i],1/sigma**2)
mu[i]<-((1-type[i])*mu_ind[individual[i]])+(type[i]*mu_ind[individual[i]]*shrink)
} 
for (j in 1:ninds){
mu_ind[j]~dunif(0,50)
}
shrink~dunif(0,3)
sigma~dnorm(1,1/4**2)}
"
writeLines(model_shrink, con="model_shrink_factor_JAGS.txt" )

parameters<-c("shrink")
adaptSteps <- 500
num_iterations<-1000
burn_in <- 1000
thinSteps<-10
num_cadenas<-3

set.seed(26)

runJagsOut <- run.jags(model="model_shrink_factor_JAGS.txt",
                       monitor=parameters,
                       data=data_list,
                       n.chains=num_cadenas ,
                       adapt=adaptSteps,
                       burnin=burn_in,
                       sample=num_iterations,
                       thin=thinSteps,
                       summarise=FALSE,
                       plots=FALSE )

write.csv(x = summary(runJagsOut),"shrink_factor_JAGS_output.csv")