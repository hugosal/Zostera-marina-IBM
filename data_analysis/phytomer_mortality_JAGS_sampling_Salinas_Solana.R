
if(!require("runjags")){
  install.packages("runjags")
  library(runjags)}

growth_2018 <- read.csv("internodes_and_branches_2018_environment.csv", header=T)
growth_full <- rbind(growth_2018,
                     read.csv("internodes_and_branches_2000_environment.csv",
                                        header=T))

# here we see the minimum  observed number of phytomers in a rhizome
min(growth_full$number_internodes_in_rhizome)

# use 2018 data for sampling
growth <- growth_2018

stages <- matrix(nrow = max(growth$number_internodes_in_rhizome), ncol = 4)
stages[, 1] <- seq(1,max(growth$number_internodes_in_rhizome), 1)
stages[, 2] <- rep(0,max(growth$number_internodes_in_rhizome))

# A rhizome with n phytomers has one phytomer of age 1, one of age 2,
# one of age 3... one of age n
# To count the number of phytomers that have p plastochrones we
# add one to the number of phytomer of age p if the number of phytomers in
# a rhizome that has n phytomers is p <= n 
for (n in (growth$number_internodes_in_rhizome)){
  row_n<-n # numer of phytomers in rhizome
  while (row_n > 0){
    stages[row_n,2] <- stages[row_n,2]+1 # add one to phytomers of age p
    row_n <- row_n-1
  }
}
colnames(stages)<-c("p","Number of observations","Not surviving","Surviving")

# calculating proportion of initial phytomer that survive to p
for (i in 1:nrow(stages)-1){
  stages[i,3]<-stages[i,2]-stages[i+1,2]#number of dying
  stages[i,4]<-stages[i,2]-stages[i,3]#number of living
}
#to calculate the proportion of phytomer that reach age p
mortality<-as.data.frame(matrix(c(stages[,1],1-(
  stages[,4]/stages[1,2])),ncol=2,byrow=F))

colnames(mortality)<-c("p","Mortality")

data_list <- list(
  age = mortality$p, 
  mortality = mortality$Mortality,
  ntotal = dim(mortality)[1])

modelo_mort <-"model {for ( i in 1:ntotal){
mortality[i] ~ dnorm(mean[i],1/sigma**2)
mean[i] <- 1/(1 + exp(beta0 + (beta1*age[i]))) }
beta0 ~ dunif(-10, 10)
beta1 ~ dunif(-10, 10)
sigma ~ dunif(0, 10)
}"

writeLines(modelo_mort, con="model_phytomer_mortality_JAGS.txt" )
parameters <- c("beta0", "beta1", "sigma")
adaptSteps <- 500
burn_in <- 5000 
num_iter <- 10000
thinSteps <- 10
num_cadenas <- 3

runJagsOut <-run.jags(model="model_phytomer_mortality_JAGS.txt",
                      monitor=parameters,
                      data=data_list,
                      n.chains=num_cadenas,
                      adapt=adaptSteps,
                      burnin=burn_in,
                      sample=num_iter,
                      thin=thinSteps,
                      summarise=FALSE,
                      plots=FALSE)

write.csv(x = summary(runJagsOut), "phytomer_mortality_JAGS_output.csv")
