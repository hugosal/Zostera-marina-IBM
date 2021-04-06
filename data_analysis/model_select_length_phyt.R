# model selection for length of new internodes submodel
if(!require("runjags")){
  install.packages("runjags")
  library(runjags)}

if(!require("rjags")){
  install.packages("rjags")
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
intercept <- rep(1, length(internode_length))

full_dataset<-data.frame(internode_length, intercept, Temp, Temp_anom, 
                          Temp2, Temp_anom_2, GHI, Hex)


variables_set <- list(c(2, 3), c(3), c(3, 5), c(3, 5, 4),
                      c(3,5,4,6), c(3, 5, 4, 6, 7), c(3, 5, 4, 6, 7, 8),
                      c(3, 5, 6, 7, 8))

dir.create("model_select_length_phyt")
for (vars in seq_along(variables_set)) {

  indep<-as.matrix(full_dataset[, variables_set[[vars]]])
  depen<-full_dataset[,1]

  data_list <- list(
    x = indep ,
    l = depen ,
    ntotal = length(depen),
    Nx = dim(indep)[2],
    ssta_squared=Temp_anom_2)

  regresion_init<-lm(depen~indep-1)
  initial <- list(
    beta = regresion_init$coef,
    betasigma = sqrt(mean(regresion_init$resid^2)))

  model_num<-"model{
  for(i in 1:ntotal){
  l[i]~dgamma(mu[i]^2/sigma[i]^2, mu[i]/sigma[i]^2)
  mu[i]<-sum(beta[1:Nx]*x[i,1:Nx])
  sigma[i]<-alfasig+betasigma*ssta_squared[i]
  }
  alfasig~dunif(0,10)
  betasigma~dnorm(0,1/2^2)
  for(j in 1:Nx){
  beta[j]~dnorm(0,1/2^2)}
  }"

  writeLines(model_num, con="model_length_selection_process.txt" )

  parameters<-c("beta", "alfasig", "betasigma" ) # which parameters to monitor
  adaptSteps <- 5000 # number of adapt iteration
  num_iter<-10000 #number or sampling iterations
  burn_in <- num_iter*0.1 #burn in
  thinSteps<-10 #thining
  num_chain<-3 #number of chains to use
  set.seed(26)
  runJagsOut <- runjags::run.jags( model="model_length_selection_process.txt",
                                   monitor=parameters,
                                   data=data_list,
                                   n.chains=num_chain,
                                   adapt=adaptSteps,
                                   burnin=burn_in,
                                   sample=num_iter,
                                   inits = initial,
                                   thin=thinSteps,
                                   summarise=FALSE,
                                   plots=FALSE,silent.jags = T) # sample
  # write each output 
  save(runJagsOut, file = paste("model_select_length_phyt/",vars, ".RData", sep = ""))
}

# next is to analyse each submodel version
variables_names <- c("depend", "Intercept", "SST", "SSTA","SST^2","SSTA^2",
                       "GHI","Hex")

model_selection <- data.frame(matrix(ncol = 5, nrow = 28))
colnames(model_selection) <- c("Model", "DIC", "Variable", "95 CI lower",
                               "95 CI upper")
r <- 1
list_files <- list.files("model_select_length_phyt/")
for ( file_number in seq_along(list_files)){
  load(paste("model_select_length_phyt/", list_files[file_number], sep = ""))
  resumen <- summary(runJagsOut)
  dic_este <- capture.output(print(extract(runJagsOut, "dic")))[9]
  dic_este <- gsub(pattern = "Penalized deviance: ", 
                   replacement = "", dic_este)
  variables_in_model <- variables_set[[file_number]]
  
  for(i in 1:nrow(resumen)){
    model_selection[r, 1] <- file_number
    model_selection[r, 2] <- dic_este
    model_selection[r, 3] <- variables_names[variables_in_model[i]]
    model_selection[r, 4] <- resumen[i, 1]
    model_selection[r, 5] <- resumen[i, 3]
    r <- r + 1
  }
}

write.csv(model_selection, file = "model_selection_table_length.csv", row.names = F)
