
if(!require("BayesFactor")){
  install.packages("BayesFactor")
  library(BayesFactor)}

# these next two files are obtained from the python code 

latin<-read.csv("latin_hipercube.csv",header = F)
out_latin<-read.csv("outputs_for_latin.csv",header = F) 

covar<-data.frame(cor(out_latin, latin, method ="pearson"))
sorted_cov<-sort(covar,decreasing = T)# compute pearson corrrelation coefficient and sort

sorted_cov<-rbind(sorted_cov,rep(0,length(sorted_cov)))

# compute BF for each correlation
for (corre in 1:length(sorted_cov)){
  sorted_cov[2,corre]<-BayesFactor::linearReg.R2stat(N=dim(latin)[1], 
                                                     p=1, R2=sorted_cov[1,corre]^2,simple = T)}

rownames(sorted_cov)<-c("Pearson correlation","Bayes Factor")

write.csv(sorted_cov,"sensitivity_table_bayes.csv")
