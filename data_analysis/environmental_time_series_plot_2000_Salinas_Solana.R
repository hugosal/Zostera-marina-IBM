
growth <- read.csv("new_internode_lenght_and_environmental.csv",header=T)
growth<-growth[complete.cases(growth), ]
dates<-unique(growth$date_stamp)
timest<-as.Date(as.POSIXct((dates - 719529)*86400, 
                           origin = "1970-01-01"))
temp<-c()
irrad<-c()
hex<-c()
anom<-c()
for (yr in dates){
  temp<-c(temp,growth$SST[growth$date_stamp==yr][1])
  irrad<-c(irrad,growth$irradiance[growth$date_stamp==yr][1])
  hex<-c(hex,growth$exposition_hours[growth$date_stamp==yr][1])
  anom<-c(anom,growth$SSTA[growth$date_stamp==yr][1])
}

postscript_file<-F

if(postscript_file){grDevices::cairo_ps("Fig_environmental_series_2000.eps", 
                                        height =5, width = 5,
                                        bg = F,fallback_resolution = 600)}else{
                                          png("Fig_environmental_series_2000.png",
                                              width = 13, height = 13, units = 'cm', res = 600) 
                                        }

plot(dates, temp, axes=F, xlab="", ylab="",type="b",lwd=2,
     col="#0096f6ff",main = "")
grid()
axis(2, ylim=c(min(temp),max(temp)),col="#0096f6ff",lwd=2,
     xlim=c(dates[1]-1,dates[32]+1),line = 0)

par(new=T)
plot(dates, irrad, axes=F, ylim=c(min(temp),max(irrad)), xlab="", ylab="", 
     type="b", main="",lwd=2,col="#8eff12ff")
axis(2, ylim=c(0,max(irrad)),lwd=2,line=2,col="#8eff12ff")

par(new=T)
plot(dates, hex, axes=F, ylim=c(min(hex),max(hex)), xlab="", ylab="", 
     type="b", main="",lwd=2,col="#e60000")
axis(4, ylim=c(min(hex),max(hex)),lwd=2,line=0,col="#e60000")

par(new=T)
plot(dates, anom, axes=F, ylim=c(min(anom),max(anom)), xlab="", ylab="", 
     type="b", main="",lwd=2,col="#b300b3")
axis(4, ylim=c(min(anom),max(anom)),lwd=3,line=2,col="#b300b3")

axis(1,dates,labels = timest,las=3)

legend(x = c(10,10),legend = c("Temperature (Â°C)",
    as.expression(bquote("Irrad (kW "~fortnight^-1~m^-2*~")")),
    "Hex","Temperature anomaly (Â°C)"), fill =c("#0096f6ff",
         "#8eff12ff","#e60000","#b300b3"),ncol = 1,inset=c(3,4),
    xpd = T)
dev.off()

write.table(data.frame(temperature=temp,irradiance=irrad,
                       temperature_anomaly=anom,exposition_hours=hex),
            "environment_data.txt")

