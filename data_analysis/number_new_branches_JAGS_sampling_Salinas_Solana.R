
if (!file.exists("number_new_branches_JAGS_output.csv")){
  
  if(!require("runjags")){
    install.packages("runjags")
    library(runjags)}
  
  growth <- read.csv("internodes_and_branches_2000.csv",header=T)
  
  n_branch<-growth$new_branches_per_branch 
  initial_branches<-growth$number_initial_branches[!(is.nan(n_branch))]
  n<-growth$new_leaves_number[!(is.nan(n_branch))]
  final_branches<-growth$number_final_branches[!(is.nan(n_branch))]
  
  data_list <- list(
    final_branches = final_branches,
    n=n,
    initial_branches=initial_branches,
    ntotal = length(initial_branches))
  
  branch_model <-"model {for (i in 1:ntotal){
  final_branches[i]~dpois(initial_branches[i]+born_branches[i]-lost_branches[i])
  born_branches[i]~dbinom(p_branch, n[i]*initial_branches[i])
  lost_branches[i]~dbinom(p_loss,initial_branches[i]-1)
  }
  p_branch~dbeta(1,1)
  p_loss~dbeta(1,1)}
  "
  writeLines(branch_model, con="model_number_new_branches_JAGS.txt")
  
  parameters<-c("p_branch", "p_loss")
  burn_in <- 5000 
  num_iters<-10000
  adaptSteps <- num_iters*0.1 
  thinSteps<-10
  num_chains<-3
  set.seed(26)
  runJagsOut <- run.jags( model="model_number_new_branches_JAGS.txt", 
                          monitor=parameters, 
                          data=data_list,  
                          n.chains=num_chains,
                          adapt=adaptSteps,
                          burnin=burn_in,
                          sample=num_iters,
                          thin=thinSteps,
                          summarise=FALSE,
                          plots=FALSE )
  
  write.csv(x = summary(runJagsOut),"number_new_branches_JAGS_output.csv")
}

