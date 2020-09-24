# ST533 MT1 Data Cleaning by HJ
# 09/21/2020

# Init
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read Data
data <- read.csv(file = "PM2.5Clean.csv", header = TRUE)

# Findout Weekends
data$weekend <- rep(FALSE, 22891)
for (ii in 1:dim(data)[1]){
  tmp <- data[ii,]
  date <- as.Date(paste(as.character(tmp$Year), as.character(tmp$Mon), as.character(tmp$Day), sep = "-"))
  mark <- weekdays(date, abbr = TRUE)
  if (mark == "Sat" | mark == "Sun"){
    data$weekend[ii] <- TRUE
  }
}

# Loop for each site
data.final <- data.frame(0, "a","Wake" ,"NC", 12, 34, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0)
names(data.final) <- c("Davg", "Site", "County", "State", "Long", "Lat", "AQI19", "AQI20", "Daqi", 
                       "W19", "W20", "E19", "E20", "dw", "de", "mV19", "mV20", "dmv", "Nobs19", "Nobs20")
lvl <- levels(data$Site.Name)
for(ii in seq_along(lvl)){
  tmp <- data[data$Site.Name == lvl[ii],]
  t19 <- tmp[tmp$Year == 2019,]
  t20 <- tmp[tmp$Year == 2020,]
  #Difference in daily PM2.5 concentration
  davg <- - mean(t19$Daily.Mean.PM2.5.Concentration) + mean(t20$Daily.Mean.PM2.5.Concentration)
  #Site's information
  site <- tmp$Site.Name[1]
  county <- tmp$COUNTY[1]
  state <- tmp$STATE[1]
  #longtitude $ Latitude
  long <- tmp$SITE_LONGITUDE[1]
  lat <- tmp$SITE_LATITUDE[1]
  #AQI
  aqi19 <- mean(t19$DAILY_AQI_VALUE)
  aqi20 <- mean(t20$DAILY_AQI_VALUE)
  daqi <- aqi19 - aqi20
  #Week vs Weekend
  w19 <- mean(t19$Daily.Mean.PM2.5.Concentration[t19$weekend == FALSE])
  w20 <- mean(t20$Daily.Mean.PM2.5.Concentration[t20$weekend == FALSE])
  dw <- - w19 + w20
  e19 <- mean(t19$Daily.Mean.PM2.5.Concentration[t19$weekend == TRUE])
  e20 <- mean(t20$Daily.Mean.PM2.5.Concentration[t20$weekend == TRUE])
  de <- - e19 + e20
  #"mean" Variance
  Nobs19 <- dim(t19)[1]
  Nobs20 <- dim(t20)[1]
  mv19 <- var(t19$Daily.Mean.PM2.5.Concentration) / Nobs19
  mv20 <- var(t20$Daily.Mean.PM2.5.Concentration) / Nobs20
  dmv <- - mv19 + mv20
  current <- data.frame(davg, site, county, state, long, lat, aqi19, aqi20, daqi, 
                        w19, w20, e19, e20, dw, de, mv19, mv20, dmv, Nobs19, Nobs20)
  names(current) <- c("Davg", "Site", "County", "State", "Long", "Lat", "AQI19", "AQI20", "Daqi", 
                         "W19", "W20", "E19", "E20", "dw", "de", "mV19", "mV20", "dmv", "Nobs19", "Nobs20")
  data.final <- rbind(data.final, current)
}
data.final <- data.final[2:dim(data.final)[1],]
data.final <- data.final[data.final$Nobs19 > 10 & data.final$Nobs20 > 10, ]
c <- cor(as.matrix(data.final[,c(1,5:20)]))

# Save the data
save(data.final, file = "envdata.RData")
