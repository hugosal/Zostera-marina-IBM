# model selection for number of new phytomers submodel

if(!require("runjags")){
  install.packages("runjags")
  library(runjags)}

if(!require("rjags")){
  install.packages("rjags")
  library(runjags)}

growth <- read.csv("internodes_and_branches_2018_environment.csv",header=T)
new_leaves<-growth$new_leaves_number 
intercept <- rep(1, length(new_leaves)) 
Temp<-growth$SST 
Temp2<-Temp^2 
Temp_anom<-growth$SSTA 
Temp_anom_2<-Temp_anom^2 
GHI<-growth$irradiance 
Hex<-growth$exposition_hours

full_dataset<-data.frame(new_leaves,intercept, Temp, Temp_anom, 
                          Temp2, Temp_anom_2, GHI, Hex)

variables_set <- list(c(2), c(2,6),  c(2,5), c(2,5,6), c(2,4), c(2, 7), c(2, 8),
                      c(2,3))

dir.create("model_select_n_phyt")
for (vars in seq_along(variables_set)) {
  
  indep<-as.matrix(full_dataset[,variables_set[[vars]]])
  depen<-full_dataset[,1]

  data_list <- list(
    x = indep ,
    new_leaves = depen ,
    ntotal = length(depen),
    Nx = dim(indep)[2])

  model_num <-"model {
  for ( i in 1:ntotal ) {
  new_leaves[i]~dpois(lambda[i])
  lambda[i]<- exp(sum(beta[1:Nx]*x[i,1:Nx]))
  }
  for(j in 1:Nx){
  beta[j]~dnorm(0,1/2^2)}
  }"

  writeLines(model_num, con="model_number_selection_process.txt" )

  parameters<-c("beta") # which parameterrs to monitor
  adaptSteps <- 5000 # number of adapt iteration
  num_iter<-10000 #number or sampling iterations
  burn_in <- num_iter*0.1 #burn in
  thinSteps<-10 #thining
  num_chain<-3 #number of chains to use
  set.seed(26)
  runJagsOut <- runjags::run.jags( model="model_number_selection_process.txt",
                          monitor=parameters,
                          data=data_list,
                          n.chains=num_chain,
                          adapt=adaptSteps,
                          burnin=burn_in,
                          sample=num_iter,
                          thin=thinSteps,
                          summarise=FALSE,
                          plots=FALSE,silent.jags = T) # sample
  # write each output 
  save(runJagsOut, file = paste("model_select_n_phyt/",vars, ".RData", sep = ""))

}

# next is to analyse each submodel version
variable_names <- c("depend", "Intercept", "SST", "SST^2","SSTA","SSTA^2",
                       "GHI","Hex")

model_selection <- data.frame(matrix(ncol = 5, nrow = 16))
colnames(model_selection) <- c("Model", "DIC", "Variable", "95 CI lower",
                               "95 CI upper")
r <- 1
list_files <- list.files("model_select_n_phyt/")
for ( fil in seq_along(list_files)){
  load(paste("model_select_n_phyt/",list_files[fil], sep = ""))
  resumen <- summary(runJagsOut)
  dic_current <- capture.output(print(extract(runJagsOut, "dic")))[9]
  dic_current<- gsub(pattern = "Penalized deviance: ", 
                             replacement = "", dic_current)
  variables_in_model <- variables_set[[fil]]
  for(i in 1:nrow(resumen)){
    model_selection[r, 1] <- fil
    model_selection[r, 2] <- dic_current
    model_selection[r, 3] <- variable_names[variables_in_model[i]]
    model_selection[r, 4] <- resumen[i, 1]
    model_selection[r, 5] <- resumen[i, 3]
    r <- r + 1
  }
}

write.csv(model_selection, file = "model_selection_table_n_phyt.csv", row.names = F)
