
if(!require("BayesFactor")){
  install.packages("BayesFactor")
  library(BayesFactor)}

want_p_values <- FALSE

# these two files are the output of the python code "sensitivity_sampling.py"

latin<-read.csv("latin_hipercube.csv",header = T)
out_latin<-read.csv("outputs_for_latin.csv",header = F) 

correl<-data.frame(cor(out_latin, latin, method ="pearson"))
bf_line <- 2
if (want_p_values) {
  bf_line <- 4
  pes <- c()
  tes <- c()
  for (var in 1:ncol(latin)){
    test <- cor.test(latin[,var],out_latin[,1])
    pes <- c(pes, round(test$p.value,3))
    tes <- c(tes, round(test$statistic,3))
  }
}

if (want_p_values) {
  correl<-rbind(correl, pes, tes)
}

sorted_corr<-correl[, colnames(sort(correl[1, ], decreasing = T))]

sorted_corr<-rbind(sorted_corr,rep(0,length(sorted_corr)))

# compute BF for each correlation


for (corre in 1:length(sorted_corr)){
  sorted_corr[bf_line, corre]<-BayesFactor::linearReg.R2stat(N=dim(latin)[1], 
                                                     p=1, 
                                                     R2=sorted_corr[1, corre]^2,
                                                     simple = T)}
if (want_p_values) {
  rownames(sorted_corr)<-c("Pearson correlation","p-value","t","Bayes Factor")
}else{
  rownames(sorted_corr)<-c("Pearson correlation", "Bayes Factor")
}

tabl<-t(sorted_corr)

write.csv(tabl,"sensitivity_table_bayes.csv")
