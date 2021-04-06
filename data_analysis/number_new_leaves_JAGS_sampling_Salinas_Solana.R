
if(!file.exists("number_new_leaves_JAGS_output.csv")){
  
  if(!require("runjags")){
    install.packages("runjags")
    library(runjags)}
    
  growth <- read.csv("internodes_and_branches_2018_environment.csv",header=T)
  new_leaves<-growth$new_leaves_number 
  
  data_list <- list(
    new_leaves = new_leaves,
    ntotal = length(new_leaves)) 
  
  model_num <-"model {
  for ( i in 1:ntotal ) {
  new_leaves[i]~dpois(lambda)
  }
  lambda~dunif(0, 10)
  }" 
  
  writeLines(model_num, con="model_number_new_leaves_JAGS.txt" )
  
  parameters<-c("lambda") # which parameterrs to monitor
  adaptSteps <- 5000 # number of adapt iteration
  num_iter<-10000 #number or sampling iterations
  burn_in <- num_iter*0.1 #burn in
  thinSteps<-10 #thining 
  num_chain<-3 #number of chains to use
  set.seed(26)
  runJagsOut <- run.jags( model="model_number_new_leaves_JAGS.txt", 
                          monitor=parameters, 
                          data=data_list,  
                          n.chains=num_chain,
                          adapt=adaptSteps,
                          burnin=burn_in, 
                          sample=num_iter,
                          thin=thinSteps,
                          summarise=FALSE,
                          plots=FALSE,silent.jags = T) # sample
  
  write.csv(x = summary(runJagsOut),"number_new_leaves_JAGS_output.csv")
}
