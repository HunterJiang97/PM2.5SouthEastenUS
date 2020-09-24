# ST533 MT1 Model Checking by HJ
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
ols <- lm(Davg ~ Long + Lat + dw ,data = data.final)
fit <- ols$fitted.values
res <- data.final$Davg-fit

df <- data.frame(long=data.final$Long,lat=data.final$Lat,Y=res)
ggplot(df, aes(long, lat)) +
  borders("state") +
  geom_point(aes(colour = Y)) +
  scale_colour_gradientn(colours = viridis(10)) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  xlab("")+ylab("")+labs(title="Observations")

s <- cbind(data.final$Long, data.final$Lat)
L     <- 30
d_max <- 10
d     <- seq(0,d_max,length=L)
d
vg <- variog(coords = s, data = res, uvec = d) 
vg$n
plot(vg) 
d <- seq(0,10,.01)
sig2 <- 0.06
tau2 <- 0.03
rho  <- 1.5
vg_fitted1 <- sig2 + tau2 - sig2*exp(-d/rho)

df2 <- data.frame(u = vg$u, v = vg$v)
df1 <- data.frame(vg = vg_fitted1, d=d)
ggplot() + 
  geom_line(aes(x = d, y = vg), data = df1, size = 1.1) + 
  geom_point(data = df2, aes(x = u, y = v), size = 2, shape = 10,  colour = "blue")

g1 <- ifelse(s[,2] < 33,1,2)
g1 <- ifelse(s[,1] > -82,g1+2,g1)
g <- g1
plot(s,col=g1)
abline(33,0,col=2)
abline(v=-82,col=3)

vg1 <- variog(coords = s[g==1,], data = res[g==1], uvec = seq(0,5,length=10)) 
vg2 <- variog(coords = s[g==2,], data = res[g==2], uvec = seq(0,5,length=10)) 
vg3 <- variog(coords = s[g==3,], data = res[g==3], uvec = seq(0,5,length=10)) 
vg4 <- variog(coords = s[g==4,], data = res[g==4], uvec = seq(0,5,length=10)) 

plot(vg1$u,vg1$v,type="l",ylim=c(0,.8))
lines(vg2$u,vg2$v,col=2)
lines(vg3$u,vg3$v,col=3)
lines(vg4$u,vg4$v,col=4)
legend("topleft",paste("Subregion",1:4),lty=1,col=1:4,ncol=2,bty="n")


df2 <- data.frame(u = vg$u, v = vg$v)
df1 <- data.frame(vg = vg_fitted1, d=d)
ggplot() + 
  geom_line(aes(x = d, y = vg), data = df1, size = 1.1) + 
  geom_point(data = df2, aes(x = u, y = v), size = 2, shape = 17) + xlab("distance") + ylab("semivariance")

region = c(rep(1, length(vg1$u)), rep(2, length(vg2$u)), rep(3, length(vg3$u)), rep(4, length(vg4$u)))
df3 <- data.frame(u = c(vg1$u,vg2$u,vg3$u,vg4$u), v = c(vg1$v,vg2$v,vg3$v,vg4$v), Region = as.factor(region))
ggplot() +
  geom_line(data = df3, aes(x = u, y = v, colour = Region, group = Region), size = 1.1) +
  ylim(0, 0.3) + xlab("distance") + ylab("semivariance")

