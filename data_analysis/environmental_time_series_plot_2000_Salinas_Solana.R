postscript_file <- T

if(postscript_file){grDevices::postscript("Fig3.eps", 
     height =5.51, width = 5.51, horizontal = F, onefile = FALSE, 
     paper = "special")
  }else{
      png("fig_environmental_time_Series.png", 
      width = 13, height = 13, units = 'cm', res = 600)
     }

colores <- rep("#666666", 4)# gray.colors(n = 4, start = 0, end = 0.4)
pchs <- c(16, 16, 16, 16)

par(oma = c(6, 2.8, 2, 0))
layout(matrix(c(1,1,5,2,2,6,3,3,7,4,4,8), nrow = 4, ncol = 3, byrow = TRUE))
second <- FALSE
for(file in c("new_internode_lenght_and_environmental.csv", 
  "new_internode_lenght_and_environmental_2018.csv" )){
 if (!second){
 par(mar = c(0.25, 2, 0.25, 0))
 }else{
 par(mar = c(0.25, 0, 0.25, 1))
 }
growth <- read.csv(file, header=T)

growth <- growth[complete.cases(growth), ]

growth <- growth[growth$length_of_internode < 50, ] # remove three outliers

dates <- unique(growth$date_stamp)

timest <- format(as.Date(as.character(as.POSIXct((dates - 719529)*86400, 
      origin = "1970-01-01")), 
   format = "%Y-%m-%d"), "%d/%m")
  
mean_lenght_new_inte <- c()
sd_lenght_new_inte <- c()
mean_number_new_phy <- c()
sd_number_new_phy <- c()
temp <- c()
irrad <- c()
hex <- c()
anom <- c()
for (yr in dates){
 monts_curren <- growth[growth$date_stamp==yr, ]
 
 mean_number_new_phy <- c(mean_number_new_phy, 
    mean(unlist(lapply(split(monts_curren$Numbrer_of_internodes_in_rhizome, 
       f = monts_curren$rhizome_number, drop = T), length))))
 
 sd_number_new_phy <- c(sd_number_new_phy, 
                       sd(unlist(lapply(split(monts_curren$Numbrer_of_internodes_in_rhizome, 
                                                f = monts_curren$rhizome_number, 
                                                drop = T), length))))
 mean_lenght_new_inte <- c(mean_lenght_new_inte, 
    mean(growth$length_of_internode[growth$date_stamp==yr], 
     na.rm = T))
 
 sd_lenght_new_inte <- c(sd_lenght_new_inte, 
                           sd(growth$length_of_internode[growth$date_stamp==yr], 
                                na.rm = T))
 
 temp <- c(temp, growth$SST[growth$date_stamp==yr][1])
 irrad <- c(irrad, growth$irradiance[growth$date_stamp==yr][1])
 hex <- c(hex, growth$exposition_hours[growth$date_stamp==yr][1])
 anom <- c(anom, growth$SSTA[growth$date_stamp==yr][1])
}

lim_t <- c(14, 27)
lim_irrad <- c(45, 250)
lim_nl <- c(1, 50)
lim_nh <- c(0.4, 5)

plot(dates, temp, axes=F, xlab="", ylab="", type="o", lwd=1, 
 col=colores[1], main = "", pch=pchs[1], ylim=lim_t)
grid()
points(dates, temp, type="o", lwd=1, col=colores[1], pch=pchs[1])
if (!second){
axis(2, lwd=1, 
 xlim=c(dates[1]-1, dates[32]+1), line = 0, 
 at = round(seq(from=lim_t[1], to=lim_t[2], length.out = 6)), cex.axis=0.8)
mtext("Water \ntemperature (°C)", side=2, line=2, cex = 0.75)
}

plot(dates, irrad, axes=F, xlab="", ylab="", 
 type="o", main="", lwd=1, pch=pchs[2], col=colores[2], 
 ylim=lim_irrad)
grid()
points(dates, irrad, type="o", lwd=1, pch=pchs[2], col=colores[2])
if (!second){
axis(2, lwd=1, line=0, 
 at = round(seq(from=lim_irrad[1], to=lim_irrad[2], length.out = 6)), 
 cex.axis=0.8)
mtext("Irradiance  (kW", side=2, line=3.5, cex = 0.75) 
mtext(as.expression(bquote(fortnight^-1~m^-2*~")")), 
      side=2, line=2, cex = 0.75) 
}

plot(dates, mean_lenght_new_inte, axes=F, ylim=lim_nl, xlab="", ylab="", 
 type="o", main="", lwd=1, pch=15, col=colores[3])
grid()
for (t in 1:length(mean_lenght_new_inte)){
  lines(x = rep(dates[t], 2), y = c(mean_lenght_new_inte[t]-sd_lenght_new_inte[t],
                                    mean_lenght_new_inte[t]+sd_lenght_new_inte[t]),
        lwd=1, col=colores[3])}
points(dates, mean_lenght_new_inte, type="o", lwd=1, pch=15, col=colores[3])

if (!second){
axis(2, at = round(seq(from=lim_nl[1], to=lim_nl[2], 
   length.out = 6)), lwd=1, line=0, cex.axis=0.8)
mtext("Mean length of \nnew internodes (mm)", side=2, line=2, cex = 0.75) 
}

plot(dates, mean_number_new_phy, axes=F, ylim=lim_nh, xlab="", ylab="", 
 type="o", main="", lwd=1, pch=15, col=colores[4])
grid()
for (t in 1:length(mean_lenght_new_inte)){
  lines(x = rep(dates[t], 2), y = c(mean_number_new_phy[t]-sd_number_new_phy[t],
                                    mean_number_new_phy[t]+sd_number_new_phy[t]),
        lwd=1, col=colores[4])}
points(dates, mean_number_new_phy, type="o", main="", lwd=1, pch=15, col=colores[4])
if (!second){
axis(2, at = round(seq(from = lim_nh[1], to =lim_nh[2], length.out = 6), digits = 1), 
 lwd=1, line=0, cex.axis=0.8)
mtext("Mean number of \n new phytomers", side=2, line=2, cex = 0.75) 
}

axis(1, dates, labels = F)
text(x = dates,  y = par("usr")[3] - 0.9,  labels = timest,  xpd = NA, 
 srt = 90,  cex = 0.7)


if (!second){
mtext("                          Date of sample", side=1, line=4)
  mtext("1999-2000", side=1, line=2.5, cex = 0.8)
}else{
  mtext("2018", side=1, line=2.5, cex = 0.8)
}
second <- !second
}

dev.off()



