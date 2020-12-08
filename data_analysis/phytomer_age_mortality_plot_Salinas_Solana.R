
growth <- read.csv("internodes_and_branches_2000.csv",header=T)

age_class<-matrix(nrow = max(growth$number_internodes_in_rhizome),ncol = 4)
age_class[,1]<-seq(1,max(growth$number_internodes_in_rhizome),1)
age_class[,2]<-rep(0,max(growth$number_internodes_in_rhizome))
for (n in (growth$number_internodes_in_rhizome)){
  ag_cls<-n
  while (ag_cls>0){
    age_class[ag_cls,2]<-age_class[ag_cls,2]+1
    ag_cls<-ag_cls-1
  }
}

colnames(age_class)<-c("Age","Number_of_observations","dead","survived")

# calculating proportion of initial phytomer that survive to p
for (i in 1:nrow(age_class)-1){
  age_class[i,3]<-age_class[i,2]-age_class[i+1,2]#number of dying
  age_class[i,4]<-age_class[i,2]-age_class[i,3]#number of living
}
#to calculate the proportion of phytomer that reach age p
cummulative_mortality<-as.data.frame(matrix(c(age_class[,1],1-(
  age_class[,4]/age_class[1,2])),ncol=2,byrow=F))
colnames(cummulative_mortality)<-c("Age","Mortality")
plot(cummulative_mortality,xlab="Age (plastochrone)",ylab="Mortality")
