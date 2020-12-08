

if(!require("BayesFactor")){
  install.packages("BayesFactor")
  library(BayesFactor)}

bayes_anov_zostera<-function(simulated_file_name,observed_file_name){
  observed <- read.csv(observed_file_name,header=T,
                       stringsAsFactors = F)
  observed<-observed[observed$Internode_length>0,]
  simulated <- read.csv(simulated_file_name,header=T,
                        stringsAsFactors = F)
  simulated<-simulated[,c(1,2,4)]
  rhizome_length <- data.frame(date=character(),
                               type=character(),
                               lengths=numeric(),
                               stringsAsFactors = F)
  dates <- intersect(unique(as.character(simulated[,1])),
                     unique(as.character(observed[,1])))# intersect dates
  dates <- dates[!is.na(dates)]
  rhizome_number <-intersect(unique(as.character(simulated[,2])),
                             unique(as.character(observed[,2])))
  rhizome_number <- rhizome_number[!is.na(rhizome_number)]
  i<-1
  for (f in dates){
    for (r in rhizome_number){
      rhizome_length[i,1]<-f
      rhizome_length[i,2]<-"Observed"
      rhizome_length[i,3]<-sum(observed[which(observed[,1]==f & observed[,2]==r),3],
                               na.rm = T)
      i<-i+1
      rhizome_length[i,1]<-f
      rhizome_length[i,2]<-"Simulated"
      rhizome_length[i,3]<-sum(simulated[which(simulated[,1]==f & simulated[,2]==r),3],
                               na.rm = T)
      i<-i+1
    }} # compute the length of each rhizome as the sum of its internodes lengths 
  
  rhizome_length<-rhizome_length[!rhizome_length$lengths==0,]
  rhizome_length<-rhizome_length[rhizome_length$lengths>1,]
  rhizome_length<-rhizome_length[!(rhizome_length[,3]==0),]
  
  rhizome_length$date<-as.factor(rhizome_length$date)
  rhizome_length$type<-as.factor(rhizome_length$type)
  
  an<-BayesFactor::anovaBF(lengths ~ date + type, data=rhizome_length,progress = F)#calculate BF
  
  result_vector<-rep(0,4) # initialize vector
  result_vector[1]<-BayesFactor::extractBF(an[1])[1]
  result_vector[2]<-BayesFactor::extractBF(an[2])[1]
  result_vector[3]<-BayesFactor::extractBF(an[3])[1]
  result_vector[4]<-BayesFactor::extractBF(an[4])[1]
  return(result_vector)
}

set.seed(26)

# for 2000
observed2000<-"observed_rhizomes_2000.csv"# the observed rhizomes file
files2000 <- list.files(path="outs2000/",full.names = T) # the directory with many simulated rhizome files

out2000<-lapply(files2000,FUN = bayes_anov_zostera,observed_file_name=observed2000)
out2000<-data.frame(matrix(unlist(out2000),nrow=length(out2000), byrow=T))

files2000 <- list.files(path="outs2000/",full.names = F) # file names 
# to make row name equal to random seed number
get_output_random_seed<-function(name){
  substr(name,start = 80,stop = 81)}

rownames(out2000)<-unlist(lapply(files2000, FUN = get_output_random_seed))

#column names
colnames(out2000)<-c("Time","Type", "Time+Type","Time+Type+TimeXType")


#for 2018
observed2018<-"observed_rhizomes_2018.csv"
files2018 <- list.files(path = "outs2018/",full.names = T)
out2018<-lapply(files2018,FUN = bayes_anov_zostera,
                observed_file_name=observed2018)
out2018<-data.frame(matrix(unlist(out2018),nrow=length(out2018), byrow = T))
files2018 <- list.files(path = "outs2018/",full.names = F)
rownames(out2018)<-unlist(lapply(files2018,FUN = get_output_random_seed))
colnames(out2018)<-c("Time","Type", "Time+Type","Time+Type+TimeXType")

# function to compute standart error
stand_error<-function(x){sd(x)/sqrt(length(x))}

# join 2000 and 2018 results in a single table, and compute the mean and SE of BF
validation_tab<-data.frame(mean=c(apply(out2000,MARGIN = 2,FUN = mean),
                                  apply(out2018,MARGIN = 2,FUN = mean)),
                           error=c(apply(out2000,MARGIN = 2,FUN = stand_error),
                                   apply(out2018,MARGIN = 2,FUN = stand_error)))
# Add a column with a specific row names
validation_tab<-cbind(data.frame(Factor=c("2000 Time",
                                          "2000 Type", "2000 Type+Type",
                                          "2000 Type+Type+TypeXTime",
                                          "2018 Time","2018 Type","2018 Type+Type",
                                          "2018 Type+Type+TypeXTime")),validation_tab)

write.csv(validation_tab,file = "validation_BF_table.csv")

