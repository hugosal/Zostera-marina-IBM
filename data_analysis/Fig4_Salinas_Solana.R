

postscript_file <- T

if(postscript_file){postscript("Fig4.eps", 
        horizontal = FALSE, onefile = FALSE, 
        paper = "special", height = 5.51, 
        width =5.51)}else{
         png("Fig4.png", 
          width = 9.5, height = 20, 
          units = 'cm', res = 600)}


par(mfrow=c(2, 1))
par(mar=c(3, 4, 2, 1))

color_palette <- c("gray0", "gray57")

data_simulated_2000 <- read.csv("outs2000/ecol_mat_nod_founding_rhizomes_2000_environment_2000_cannal_200m_broad_4m_prof_31.csv", 
        header=T, stringsAsFactors = F)
data_simulated_2000 <- data_simulated_2000[, c(1, 2, 4)] 
data_observed_2000 <- read.csv("observed_rhizomes_2000.csv", header=T, stringsAsFactors = F)
data_observed_2000 <- data_observed_2000[data_observed_2000$Internode_length>0, ]

rhizome_length <- data.frame(date=character(), tipo=character(), 
        longitud=numeric(), stringsAsFactors = F)
dates <- intersect(unique(as.character(data_simulated_2000[, 1])), 
     unique(as.character(data_observed_2000[, 1])))
dates <- dates[!is.na(dates)]
rhizomes_numbers <- intersect(unique(as.character(data_simulated_2000[, 2])), 
        unique(as.character(data_observed_2000[, 2])))
rhizomes_numbers <- rhizomes_numbers[!is.na(rhizomes_numbers)]

i <- 1
for (f in dates){
 for (r in rhizomes_numbers){
 rhizome_length[i, 1] <- f
 rhizome_length[i, 2] <- "observed"
 rhizome_length[i, 3] <- sum(data_observed_2000[which(data_observed_2000[, 1]==f & 
              data_observed_2000[, 2]==r), 3], 
        na.rm = T)
 i <- i+1
 rhizome_length[i, 1] <- f
 rhizome_length[i, 2] <- "simulated"
 rhizome_length[i, 3] <- sum(data_simulated_2000[which(data_simulated_2000[, 1]==f &
               data_simulated_2000[, 2]==r), 3], 
        na.rm = T)
 i <- i+1
 }}

rhizome_length <- rhizome_length[!rhizome_length$longitud==0, ]
rhizome_length <- rhizome_length[rhizome_length$longitud>1, ]
rhizome_length <- rhizome_length[!(rhizome_length[, 3]==0), ]
rhizome_length$date <- as.factor(rhizome_length$date)
rhizome_length$tipo <- as.factor(rhizome_length$tipo)

mean_length <- data.frame()
i <- 1
for (d in dates){
 mean_length[i, 1] <- as.character(d)
 mean_length[i, 2] <- mean(rhizome_length$longitud[
 which(rhizome_length$date==d &
   rhizome_length$tipo=="observed")])
 mean_length[i, 3] <- mean(rhizome_length$longitud[
 which(rhizome_length$date==d &
   rhizome_length$tipo=="simulated")]) 

 mean_length[i, 4] <- sd(rhizome_length$longitud[
 which(rhizome_length$date==d &
   rhizome_length$tipo=="observed")]) 
 mean_length[i, 5] <- sd(rhizome_length$longitud[
 which(rhizome_length$date==d &
   rhizome_length$tipo=="simulated")])
 i <- i+1
}
colnames(mean_length) <- c("date", "observed", "simulated", "sd_obs", "sd_sim")
tiempos <- length(dates)

plot(1:tiempos, rep(1, tiempos), 
  type='n', xaxt = "n", ylim=c(50, 420), 
  ylab="Rhizome length (mm)", xlab="", main="Year 2000")
grid()
for (t in 1:tiempos){
 obs_sup <- mean_length$observed[t]+mean_length$sd_obs[t]
 obs_inf <- mean_length$observed[t]-mean_length$sd_obs[t]
 lines(x = rep(t, 2)-0.05, y = c(obs_inf, obs_sup), col=color_palette[1], lty=1)
 
 sim_sup <- mean_length$simulated[t]+mean_length$sd_sim[t]
 sim_inf <- mean_length$simulated[t]-mean_length$sd_sim[t]
 lines(x = rep(t, 2)+0.05, y = c(sim_inf, sim_sup), col=color_palette[2], lty=1)
 
}

lines(1:tiempos-0.05, mean_length$observed, 
  col=color_palette[1], type = "l")
lines(1:tiempos+0.05, mean_length$simulated, 
  col=color_palette[2], type = "l")

points(1:tiempos-0.05, mean_length$observed, 
  col=color_palette[1], pch=15, type = "b")
points(1:tiempos+0.05, mean_length$simulated, 
  col=color_palette[2], pch=17, type = "b")

axis(1, at = 1:length(dates), labels =unique(substr(rhizome_length$date, 1, 5)), 
  las = 2, cex.axis=0.6)
legend(x="top", ncol = 2, legend=c("Observed", "Simulated"), 
  col = c(color_palette[1], color_palette[2]), cex=1, pch=c(15, 17), lty = 1)


data_simulated_2018 <- read.csv("outs2018/ecol_mat_nod_founding_rhizomes_2018_environment_2018_cannal_200m_broad_4m_prof_32.csv", header=T, stringsAsFactors = F)
data_simulated_2018 <- data_simulated_2018[, c(1, 2, 4)] 
data_observed_2018 <- read.csv("observed_rhizomes_2018.csv", header=T, stringsAsFactors = F)
data_observed_2018 <- data_observed_2018[data_observed_2018$Internode_length>0, ]
rhizome_length_2018 <- data.frame(date=character(), tipo=character(), longitud=numeric(), 
         stringsAsFactors = F)
dates_2018 <- intersect(unique(as.character(data_simulated_2018[, 1])), 
      unique(as.character(data_observed_2018[, 1])))
dates_2018 <- dates_2018[!is.na(dates_2018)]
rhizomes_numbers <- intersect(unique(as.character(data_simulated_2018[, 2])), 
        unique(as.character(data_observed_2018[, 2])))
rhizomes_numbers <- rhizomes_numbers[!is.na(rhizomes_numbers)]

i <- 1
for (d in dates_2018){
 for (r in rhizomes_numbers){
 rhizome_length_2018[i, 1] <- d
 rhizome_length_2018[i, 2] <- "observed"
 rhizome_length_2018[i, 3] <- sum(data_observed_2018[which(data_observed_2018[, 1]==d & 
                data_observed_2018[, 2]==r), 3], 
         na.rm = T)
 i <- i+1
 rhizome_length_2018[i, 1] <- d
 rhizome_length_2018[i, 2] <- "simulated"
 rhizome_length_2018[i, 3] <- sum(data_simulated_2018[which(data_simulated_2018[, 1]==d & 
                data_simulated_2018[, 2]==r), 3], 
         na.rm = T)
 i <- i+1
 }}


rhizome_length_2018 <- rhizome_length_2018[!rhizome_length_2018$longitud==0, ]

rhizome_length_2018 <- rhizome_length_2018[rhizome_length_2018$longitud>1, ]

rhizome_length_2018 <- rhizome_length_2018[!(rhizome_length_2018[, 3]==0), ]


rhizome_length_2018$date <- as.factor(rhizome_length_2018$date)
rhizome_length_2018$tipo <- as.factor(rhizome_length_2018$tipo)

mean_length <- data.frame()
i <- 1
for (f in dates_2018){
 mean_length[i, 1] <- as.character(f)
 mean_length[i, 2] <- mean(rhizome_length_2018$longitud[
 which(rhizome_length_2018$date==f &
   rhizome_length_2018$tipo=="observed")]) 
 mean_length[i, 3] <- mean(rhizome_length_2018$longitud[
 which(rhizome_length_2018$date==f &
   rhizome_length_2018$tipo=="simulated")]) 

 mean_length[i, 4] <- sd(rhizome_length_2018$longitud[
 which(rhizome_length_2018$date==f &
   rhizome_length_2018$tipo=="observed")])
 mean_length[i, 5] <- sd(rhizome_length_2018$longitud[
 which(rhizome_length_2018$date==f &
   rhizome_length_2018$tipo=="simulated")])
 i <- i+1
}
colnames(mean_length) <- c("date", "observed", "simulated", "sd_obs", "sd_sim")
times2018 <- length(dates_2018)


plot(1:times2018, rep(1, times2018), 
  type='n', xaxt = "n", ylim=c(50, 900), 
  ylab="Rhizome length (mm)", xlab="", main="Year 2018" )
grid()

for (t in 1:times2018){
 obs_sup <- mean_length$observed[t]+mean_length$sd_obs[t]
 obs_inf <- mean_length$observed[t]-mean_length$sd_obs[t]
 lines(x = rep(t, 2)-0.02, y = c(obs_inf, obs_sup), col=color_palette[1], lty=1)

 sim_sup <- mean_length$simulated[t]+mean_length$sd_sim[t]
 sim_inf <- mean_length$simulated[t]-mean_length$sd_sim[t]
 lines(x = rep(t, 2)+0.02, y = c(sim_inf, sim_sup), col=color_palette[2], lty=1)
 
}

lines(1:times2018-0.02, mean_length$observed, 
  col=color_palette[1], type = "l")
lines(1:times2018+0.02, mean_length$simulated, 
  col=color_palette[2], type = "l")

points(1:times2018-0.02, mean_length$observed, 
  col=color_palette[1], pch=15, type = "b")
points(1:times2018+0.02, mean_length$simulated, 
  col=color_palette[2], pch=17, type = "b")

axis(1, at = 1:length(dates), labels=unique(substr(rhizome_length$date, 1, 5)), 
  las = 2, cex.axis=0.6)

legend(x="top", ncol = 2, legend=c("Observed", "Simulated"), 
  col = c(color_palette[1], color_palette[2]), cex=1, pch=c(15, 17), lty=1)
dev.off()
