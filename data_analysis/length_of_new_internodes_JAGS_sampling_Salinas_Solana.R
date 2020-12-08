
if(!require("runjags")){
  install.packages("runjags")
  library(runjags)}

growth <- read.csv("new_internode_lenght_and_environmental.csv",header=T)
growth<-growth[complete.cases(growth),]

zlong<-growth$length_of_internode
GHI<-growth$irradiance[zlong>0 & zlong<50]
Hex<-growth$exposition_hours[zlong>0 & zlong<50]
Temp<-growth$SST[zlong>0 & zlong<50]
Temp_anom<-growth$SSTA[zlong>0 & zlong<50]
internode_length<-zlong[zlong>0 & zlong<50]
Temp2<-Temp^2
Temp_anom_2<-Temp_anom^2

my_data_frame<-data.frame(internode_length,Temp,Temp2,Temp_anom_2,GHI,Hex)
indep<-as.matrix(my_data_frame[,2:length(my_data_frame)],
                 ncol=dim(my_data_frame)[2])
depen<-my_data_frame[,1]

data_list <- list(
  x = indep ,
  l = depen ,
  Nx = dim(indep)[2] ,
  ntotal = dim(indep)[1],
  ssta_squared=Temp_anom_2)

model_length <-"model{
for(i in 1:ntotal){
l[i]~dgamma(mu[i]^2/sigma[i]^2,mu[i]/sigma[i]^2)
mu[i]<-sum(beta[1:Nx]*x[i,1:Nx])
sigma[i]<-alfasig+betasigma*ssta_squared[i]
}
alfasig~dunif(0,10)
betasigma~dnorm(0,1/2^2)
for(j in 1:Nx){
beta[j]~dnorm(0,1/2^2)}
}"

writeLines(model_length , con="model_length_new_internodes_JAGS.txt" )

parameters <- c("beta","betasigma","alfasig") 
burn_in <- 5000 
num_iter<-5000
adaptSteps <-num_iter*0.1
thinSteps<-10
num_chains<-3

regresion_init<-lm(depen~indep-1)
initial <- list(
  beta = regresion_init$coef,        
  betasigma = sqrt(mean(regresion_init$resid^2)))

set.seed(26)

runJagsOut <- run.jags(model="model_length_new_internodes_JAGS.txt" ,
                       monitor=parameters ,
                       data=data_list ,
                       n.chains=num_chains ,
                       adapt=adaptSteps ,
                       burnin=burn_in , 
                       sample=num_iter ,
                       thin=thinSteps ,inits = initial,
                       summarise=FALSE ,
                       plots=FALSE )

#write.csv(x = summary(runJagsOut),"model_length_new_internodes_JAGS_output.csv")
