# ST533 MT1 Data Desc by HJ
# 09/22/2020

# Init
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(maps)
library(viridis)

# Read Data
load("envdata.RData")
c <- cor(as.matrix(data.final[,c(1,5:20)]))

# State
p <- ggplot(data.final, aes(x = Davg, fill = State)) +
      geom_density(alpha = 0.3)
p

# Weekend ?
df <- data.frame(Average = c(data.final$W19, data.final$E19, data.final$W20, data.final$E20), 
                 Type = rep(c(rep("Week day", 120), rep("Weekend", 120)),2),
                 Year = c(rep("2019", 240), rep("2020", 240)))
p1 <- ggplot(df, aes(x = Average, fill = Type, group = Type)) + 
  geom_density(alpha = 0.3) +
  facet_grid(rows = vars(Year))
p1

# PM2.5
ggplot(data.final, aes(Long, Lat)) +
  borders("state") +
  geom_point(aes(colour = Davg)) +
  scale_colour_gradientn(colours = viridis(10)) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45,xlim=c(-87.5,-75),ylim=c(25,40)) +
  xlab("")+ylab("")
