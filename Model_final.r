# ST533 MT1 Final Model by HJ
# 09/22/2020

# Init
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1)

# Librarys
library(ggplot2)
library(geoR)
library(knitr)
library(maps)
library(viridis)

# Read Data
load("envdata.RData")
Y <- data.final$Davg
n <- length(Y)
lon <- data.final$Long
lat <- data.final$Lat
s <- cbind(lon, lat)

# FINAL MODEL
X <- cbind(lon, lat, data.final$dw)
  
fit_mle  <- likfit(data=Y,trend= ~X,coords=s,
                    fix.nugget=FALSE,nugget=0.1,
                    cov.model="exponential",
                    ini = c(0.1, 2),messages=FALSE)
summary(fit_mle)

# DATA @dw0
s01 <- seq(-135,-65,0.1)
s02 <- seq(25,50,0.1)
s0  <- as.matrix(expand.grid(s01,s02))
inES <- map.where("state", s0[,1], s0[,2])
ind <- which(inES == "florida" | inES == "georgia" | inES == "virginia:chesapeake" |
              inES == "virginia:main" | inES == "north carolina:knotts" | inES == "north carolina:main" | inES == "south carolina")
s0 <- s0[ind,]
dw0 <- rep(NA, length(ind))
dw <- data.final$dw
for (ii in 1:length(ind)){
  dis2 <- rep(NA, 120)
  for (jj in 1:120){
    dis2[jj] <- (s[jj,1] - s0[ii,1])^2 + (s[jj,2] - s0[ii,2])^2
  }
  loc <- sort(dis2, method = "sh", index = TRUE)$ix[1:5]
  weights <- 1 / dis2[loc]
  weights <- weights / sum(weights)
  dw0[ii] <- sum(dw[loc] * weights)
}
df <- data.frame(dw0, lon = s0[,1], lat = s0[,2])
ggplot(df, aes(lon, lat)) +
  borders("state") +
  geom_point(aes(colour = dw0)) +
  scale_colour_gradientn(colours = viridis(10)) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45,xlim=c(-87.5,-75),ylim=c(25,40)) +
  xlab("")+ylab("")+labs(title="Difference in Weekday Concentration (PM2.5 2019-2020)")
X0 <- cbind(s0, dw0)


# Prediction
tick <- proc.time()
pred <- krige.conv(data=Y,coords=s, # Describe training data
                   locations=s0,    # Describe prediction sites
                   krige=krige.control(trend.d = ~X,  # Covariates at s
                                       trend.l = ~X0, # Covariates at s0
                                       cov.model="exponential",
                                       beta=fit_mle$beta,
                                       cov.pars=fit_mle$cov.pars,
                                       nugget=fit_mle$nugget))
tock <- proc.time()
tock-tick # time in seconds
Yhat <- pred$predict         # Kriging predictions
se   <- sqrt(pred$krige.var) # Kriging prediction standard deviations
PI   <- cbind(Yhat-1.96*se,Yhat+1.96*se) # Prediction intervals
PI[1:5,]

df <- data.frame(long=s0[,1],lat=s0[,2],Y=Yhat)
ggplot(df, aes(long, lat)) +
  borders("state") +
  geom_raster(aes(fill = Y)) +
  scale_fill_gradientn(colours = viridis(10))+
  xlab("")+ylab("")+labs(title="Predicted ozone (ppm), 2014")+
  coord_fixed(xlim=c(-87.5,-75),ylim=c(25,40))

# DATA @dw1
s1 <- s0
dw1 <- rep(NA, length(ind))
dX <- cbind(lon, lat, lon*lon, lat*lat, lon*lat)
lon1 <- s1[,1]
lat1 <- s1[,2]
X1 <- cbind(lon1, lat1, lon1*lon1, lat1*lat1, lon1*lat1)

dw1_mle  <- likfit(data=dw,trend= ~dX,coords=s,
                   fix.nugget=FALSE,nugget=0.1,
                   cov.model="exponential",
                   ini = c(0.1, 2),messages=FALSE)
dew_pred <- krige.conv(data=dw,coords=s, # Describe training data
                   locations=s1,    # Describe prediction sites
                   krige=krige.control(trend.d = ~dX,  # Covariates at s
                                       trend.l = ~X1, # Covariates at s0
                                       cov.model="exponential",
                                       beta=dw1_mle$beta,
                                       cov.pars=dw1_mle$cov.pars,
                                       nugget=dw1_mle$nugget))
dw1 <- dew_pred$predict
se1 <- sqrt(dew_pred$krige.var)
df <- data.frame(se1 = se1, lon = s0[,1], lat = s0[,2])
ggplot(df, aes(lon, lat)) +
  borders("state") +
  geom_point(aes(colour = se1)) +
  scale_colour_gradientn(colours = viridis(10)) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45,xlim=c(-87.5,-75),ylim=c(25,40)) +
  xlab("")+ylab("")

X1 <- cbind(s1, dw1)

# Prediction
tick <- proc.time()
pred <- krige.conv(data=Y,coords=s, # Describe training data
                   locations=s1,    # Describe prediction sites
                   krige=krige.control(trend.d = ~X,  # Covariates at s
                                       trend.l = ~X1, # Covariates at s0
                                       cov.model="exponential",
                                       beta=fit_mle$beta,
                                       cov.pars=fit_mle$cov.pars,
                                       nugget=fit_mle$nugget))
tock <- proc.time()
tock-tick # time in seconds
Yhat <- pred$predict         # Kriging predictions
se   <- sqrt(pred$krige.var) # Kriging prediction standard deviations
PI   <- cbind(Yhat-1.96*se,Yhat+1.96*se) # Prediction intervals
PI[1:5,]

df <- data.frame(long=s0[,1],lat=s0[,2],Y=Yhat)
ggplot(df, aes(long, lat)) +
  borders("state") +
  geom_raster(aes(fill = Y)) +
  scale_fill_gradientn(colours = viridis(10))+
  xlab("")+ylab("")+
  coord_fixed(xlim=c(-87.5,-75),ylim=c(25,40))

df <- data.frame(long=s0[,1],lat=s0[,2],SE=se)
ggplot(df, aes(long, lat)) +
  borders("state") +
  geom_raster(aes(fill = SE)) +
  scale_fill_gradientn(colours = viridis(10))+
  xlab("")+ylab("")+
  coord_fixed(xlim=c(-87.5,-75),ylim=c(25,40))

Yhat1 <- Yhat - 1.96*se
Yhat2 <- Yhat + 1.96*se
Yext <- rep(0, length(Yhat1))
Yext[Yhat1 > 0] <- 1
Yext[Yhat2 < 0] <- -1
df <- data.frame(long=s0[,1],lat=s0[,2],Extreme=Yext)
ggplot(df, aes(long, lat)) +
  borders("state") +
  geom_raster(aes(fill = Extreme)) +
  scale_fill_gradientn(colours = viridis(10))+
  xlab("")+ylab("")+
  coord_fixed(xlim=c(-87.5,-75),ylim=c(25,40))

Yext <- Yhat
lim <- quantile(Yext, c(0.1, 0.9))
lim
Yext[Yext < lim[1]] = -30
Yext[Yext > lim[1] & Yext < lim[2]] = -29
Yext[Yext > lim[2]] = -28
Yext[Yext == -30] = 1
Yext[Yext == -29] = 2
Yext[Yext == -28] = 3
df <- data.frame(long=s0[,1],lat=s0[,2],Extreme=Yext)
ggplot(df, aes(long, lat)) +
  borders("state") +
  geom_raster(aes(fill = Extreme)) +
  scale_fill_gradientn(colours = viridis(10))+
  xlab("")+ylab("")+
  coord_fixed(xlim=c(-87.5,-75),ylim=c(25,40))
