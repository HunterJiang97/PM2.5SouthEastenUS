# ST533 MT1 Model Comparison by HJ
# 09/22/2020

# Init
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1)

# Librarys
library(ggplot2)
library(geoR)
library(knitr)

# Read Data
load("envdata.RData")
Y <- data.final$Davg
n <- length(Y)
lon <- data.final$Long
lat <- data.final$Lat
s <- cbind(lon, lat)

# 3-fold CV
K <- 5
fold   <- sample(1:K,n,replace=TRUE)
Results <- matrix(NA, ncol = 3, nrow = 12)

# MODEL 1
X <- cbind(lon, lat)
Yhat1 <- Yhat2 <- Yhat3 <- Yhat4 <- rep(NA,n)
for(k in 1:K){
  train     <- fold != k
  test      <- fold == k
  
  fit_mle1  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=FALSE,nugget=0.1,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)
  fit_mle2  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=TRUE,nugget=0,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)  
  
  pred1 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle1$beta,
                                          cov.pars=fit_mle1$cov.pars,
                                          nugget=fit_mle1$nugget))
  
  pred2 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle2$beta,
                                          cov.pars=fit_mle2$cov.pars,
                                          nugget=fit_mle2$nugget))
  Yhat1[test] <- pred1$predict
  Yhat2[test] <- pred2$predict
}   

MSE1   <- mean((Y-Yhat1)^2)
MSE2   <- mean((Y-Yhat2)^2)

MAE1   <- mean(abs(Y-Yhat1))
MAE2   <- mean(abs(Y-Yhat2))


COR1   <- cor(Y,Yhat1)
COR2   <- cor(Y,Yhat2)


out <- rbind(c(MSE1,MSE2),
             c(MAE1,MAE2),
             c(COR1,COR2))

colnames(out) <- c("expn","exp")
rownames(out) <- c("MSE","MAE", "COR")

Results[1:2,] <- t(out)

# MODEL 2
X <- cbind(lon, lat, data.final$dw)
Yhat1 <- Yhat2 <- Yhat3 <- Yhat4 <- rep(NA,n)
for(k in 1:K){
  train     <- fold != k
  test      <- fold == k
  
  fit_mle1  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=FALSE,nugget=0.1,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)
  fit_mle2  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=TRUE,nugget=0,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)  
  
  pred1 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle1$beta,
                                          cov.pars=fit_mle1$cov.pars,
                                          nugget=fit_mle1$nugget))
  
  pred2 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle2$beta,
                                          cov.pars=fit_mle2$cov.pars,
                                          nugget=fit_mle2$nugget))
  
  Yhat1[test] <- pred1$predict
  Yhat2[test] <- pred2$predict
}   

MSE1   <- mean((Y-Yhat1)^2)
MSE2   <- mean((Y-Yhat2)^2)

MAE1   <- mean(abs(Y-Yhat1))
MAE2   <- mean(abs(Y-Yhat2))


COR1   <- cor(Y,Yhat1)
COR2   <- cor(Y,Yhat2)


out <- rbind(c(MSE1,MSE2),
             c(MAE1,MAE2),
             c(COR1,COR2))

colnames(out) <- c("expn","exp")
rownames(out) <- c("MSE","MAE", "COR")

Results[3:4,] <- t(out)


# MODEL 2
X <- cbind(lon, lat, data.final$de)
Yhat1 <- Yhat2 <- Yhat3 <- Yhat4 <- rep(NA,n)
for(k in 1:K){
  train     <- fold != k
  test      <- fold == k
  
  fit_mle1  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=FALSE,nugget=0.1,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)
  fit_mle2  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=TRUE,nugget=0,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)  
  
  pred1 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle1$beta,
                                          cov.pars=fit_mle1$cov.pars,
                                          nugget=fit_mle1$nugget))
  
  pred2 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle2$beta,
                                          cov.pars=fit_mle2$cov.pars,
                                          nugget=fit_mle2$nugget))
  
  Yhat1[test] <- pred1$predict
  Yhat2[test] <- pred2$predict
}   

MSE1   <- mean((Y-Yhat1)^2)
MSE2   <- mean((Y-Yhat2)^2)

MAE1   <- mean(abs(Y-Yhat1))
MAE2   <- mean(abs(Y-Yhat2))


COR1   <- cor(Y,Yhat1)
COR2   <- cor(Y,Yhat2)


out <- rbind(c(MSE1,MSE2),
             c(MAE1,MAE2),
             c(COR1,COR2))

colnames(out) <- c("expn","exp")
rownames(out) <- c("MSE","MAE", "COR")

Results[5:6,] <- t(out)

# MODEL 4
isFL <- isGA <- isSC <- isVA <- rep(0, n)
isFL[data.final$State == "Florida"] <- 1
isGA[data.final$State == "Georgia"] <- 1
isSC[data.final$State == "South Carolina"] <- 1
isVA[data.final$State == "Virginia"] <- 1
X <- cbind(lon, lat, isFL, isGA, isSC, isVA)
Yhat1 <- Yhat2 <- Yhat3 <- Yhat4 <- rep(NA,n)
for(k in 1:K){
  train     <- fold != k
  test      <- fold == k
  
  fit_mle1  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=FALSE,nugget=0.1,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)
  fit_mle2  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=TRUE,nugget=0,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)  
  
  pred1 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle1$beta,
                                          cov.pars=fit_mle1$cov.pars,
                                          nugget=fit_mle1$nugget))
  
  pred2 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle2$beta,
                                          cov.pars=fit_mle2$cov.pars,
                                          nugget=fit_mle2$nugget))
  
  Yhat1[test] <- pred1$predict
  Yhat2[test] <- pred2$predict
}   

MSE1   <- mean((Y-Yhat1)^2)
MSE2   <- mean((Y-Yhat2)^2)

MAE1   <- mean(abs(Y-Yhat1))
MAE2   <- mean(abs(Y-Yhat2))


COR1   <- cor(Y,Yhat1)
COR2   <- cor(Y,Yhat2)


out <- rbind(c(MSE1,MSE2),
             c(MAE1,MAE2),
             c(COR1,COR2))

colnames(out) <- c("expn","exp")
rownames(out) <- c("MSE","MAE", "COR")

Results[7:8,] <- t(out)

# MODEL 5
isFL <- isGA <- isSC <- isVA <- rep(0, n)
isFL[data.final$State == "Florida"] <- 1
isGA[data.final$State == "Georgia"] <- 1
isSC[data.final$State == "South Carolina"] <- 1
isVA[data.final$State == "Virginia"] <- 1
X <- cbind(lon, lat, data.final$dw, isFL, isGA, isSC, isVA)
Yhat1 <- Yhat2 <- Yhat3 <- Yhat4 <- rep(NA,n)
for(k in 1:K){
  train     <- fold != k
  test      <- fold == k
  
  fit_mle1  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=FALSE,nugget=0.1,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)
  fit_mle2  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=TRUE,nugget=0,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)  
  
  pred1 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle1$beta,
                                          cov.pars=fit_mle1$cov.pars,
                                          nugget=fit_mle1$nugget))
  
  pred2 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle2$beta,
                                          cov.pars=fit_mle2$cov.pars,
                                          nugget=fit_mle2$nugget))
  
  Yhat1[test] <- pred1$predict
  Yhat2[test] <- pred2$predict
}   

MSE1   <- mean((Y-Yhat1)^2)
MSE2   <- mean((Y-Yhat2)^2)

MAE1   <- mean(abs(Y-Yhat1))
MAE2   <- mean(abs(Y-Yhat2))


COR1   <- cor(Y,Yhat1)
COR2   <- cor(Y,Yhat2)


out <- rbind(c(MSE1,MSE2),
             c(MAE1,MAE2),
             c(COR1,COR2))

colnames(out) <- c("expn","exp")
rownames(out) <- c("MSE","MAE", "COR")

Results[9:10,] <- t(out)

# MODEL 6
isFL <- isGA <- isSC <- isVA <- rep(0, n)
isFL[data.final$State == "Florida"] <- 1
isGA[data.final$State == "Georgia"] <- 1
isSC[data.final$State == "South Carolina"] <- 1
isVA[data.final$State == "Virginia"] <- 1
X <- cbind(lon, lat, data.final$de, isFL, isGA, isSC, isVA)
Yhat1 <- Yhat2 <- Yhat3 <- Yhat4 <- rep(NA,n)
for(k in 1:K){
  train     <- fold != k
  test      <- fold == k
  
  fit_mle1  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=FALSE,nugget=0.1,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)
  fit_mle2  <- likfit(data=Y[train],trend= ~X[train,],coords=s[train,],
                      fix.nugget=TRUE,nugget=0,
                      cov.model="exponential",
                      ini = c(0.1, 2),messages=FALSE)  
  
  pred1 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle1$beta,
                                          cov.pars=fit_mle1$cov.pars,
                                          nugget=fit_mle1$nugget))
  
  pred2 <- krige.conv(data=Y[train],coords=s[train,], # Describe training data
                      locations=s[test,],              # Describe prediction sites
                      krige=krige.control(trend.d = ~X[train,],  # Covariates at s
                                          trend.l = ~X[test,],   # Covariates at s0
                                          cov.model="exponential",
                                          beta=fit_mle2$beta,
                                          cov.pars=fit_mle2$cov.pars,
                                          nugget=fit_mle2$nugget))
  
  Yhat1[test] <- pred1$predict
  Yhat2[test] <- pred2$predict
}   

MSE1   <- mean((Y-Yhat1)^2)
MSE2   <- mean((Y-Yhat2)^2)

MAE1   <- mean(abs(Y-Yhat1))
MAE2   <- mean(abs(Y-Yhat2))


COR1   <- cor(Y,Yhat1)
COR2   <- cor(Y,Yhat2)


out <- rbind(c(MSE1,MSE2),
             c(MAE1,MAE2),
             c(COR1,COR2))

colnames(out) <- c("expn","exp")
rownames(out) <- c("MSE","MAE", "COR")

Results[11:12,] <- t(out)

round(Results,5)

