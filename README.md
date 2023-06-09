# iCF (iterative Causal Forest)

**Citation**

Wang T, Keil AP, Kim S, Wyss R, Htoo PT, Funk MJ, Buse JB, Kosorok MR, Stürmer T. Iterative Causal Forest: A Novel Algorithm for Subgroup Identification. Am J Epidemiol. 2023 (In Press).

**R packages recommended**
```{r packages, include=FALSE}
MASS
grf
tidyverse
rlang
rlist
plyr
caret
caTools
listdtr
randomForest
ggplot2
ggridges
data.table
grid
broom
rstatix
DMwR
knitr
Rfast
spaMM
```
**installation**
```{}
install.packages('devtools')
devtools::install_github("tianshengwang/iCF")
library(iCF)
```
**Data Simulation**
```{}
  nstudy = 10000
  b0 <- 0
  b1 <- 0.8
  b2 <- -0.25
  b3 <- 0.6
  b4 <- -0.4
  b5 <- -0.8
  b6 <- -0.5
  b7 <- 0.7
  a0 <- -3.85
  a1 <- 0.3
  a2 <- -0.36
  a3 <- -0.73
  a4 <- -0.2
  a5 <- 0.71
  a6 <- -0.19
  a7 <- 0.26
  g1 <- -0.4 
  
  F.sample.cor <- function(X, rho) {
    Y <- (rho * (X - mean(X)))/sqrt(var(X)) + sqrt(1 - rho^2) * rnorm(length(X))
    return(Y)
  }

  X1 <- rnorm(nstudy, mean=0, sd=1)
  X2 <- rnorm(nstudy, mean=0, sd=1)
  X3 <- rnorm(nstudy, mean=0, sd=1)
  X4 <- rnorm(nstudy, mean=0, sd=1)
  X5 <- F.sample.cor(X1, 0.2)
  X6 <- F.sample.cor(X2, 0.9)
  X7 <- rnorm(nstudy, mean=0, sd=1)
  X8 <- F.sample.cor(X3, 0.2)
  X9 <- F.sample.cor(X4, 0.9)
  X10 <- rnorm(nstudy, mean=0, sd=1)

  X1 <- ifelse(X1 > mean(X1), 1, 0)
  X3 <- ifelse(X3 > mean(X3), 1, 0)
  X5 <- ifelse(X5 > mean(X5), 1, 0)
  X6 <- ifelse(X6 > mean(X6), 1, 0)
  X8 <- ifelse(X8 > mean(X8), 1, 0)
  X9 <- ifelse(X9 > mean(X9), 1, 0)
  
  PS <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7) ))^-1 #true propensity score
  
  W <- rbinom(nstudy,1,z.a_trueps) #treatment assignment
  
  Y = a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*a + rnorm(nstudy,0,1) + 0.4*a*X3 
  
  dat <<- as.data.frame(cbind(W, Y, X1, X2, X3 ,X4, X5, X6, X7, X8, X9, X10))
  
``` 
**Prepare data**
```{}
 vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )
 intTRUE <- "Unknown"
 X <- dat[,vars_forest]
 Y <- as.vector( as.numeric( dat[,"Y"] ) )
 W <- as.vector( as.numeric( dat[,"W"] ) )
 
 cf_raw_key.tr <- CF_RAW_key(Train, 1, "non-hd", hdpct=0.90) 
 Y.hat  <<- cf_raw_key.tr$Y.hat                 
 W.hat  <<- cf_raw_key.tr$W.hat                 
 HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    
 varimp_cf  <- cf_raw_key.tr$varimp_cf          
 length(W.hat); 
 selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx 

 PlotVI(varimp_cf, "Variable importance")
 GG_VI(varimp_cf, "Variable importance" )
 
 vars_catover2 <<- NA
 
 P_threshold <<- 0.1
  
```

**Run iCF**
```{}
D2_MLS=MinLeafSizeTune(denominator=50, treeNo = 1000, iterationNo=100, "D2")
D2_MLS$depth_mean
D2_MLS$depth_gg
![alt text](/image/D2 MLS tune.png)

D3_MLS=MinLeafSizeTune(denominator=80, treeNo = 1000, iterationNo=100, "D3")
D3_MLS$depth_mean
D3_MLS$depth_gg

D4_MLS=MinLeafSizeTune(denominator=110, treeNo = 1000, iterationNo=100, "D4")
D4_MLS$depth_mean
D4_MLS$depth_gg

D5_MLS=MinLeafSizeTune(denominator=155, treeNo = 1000, iterationNo=100, "D5")
D5_MLS$depth_mean
D5_MLS$depth_gg

leafsize <<- list(D5=155, D4=110, D3=80, D2=50)

iCFCV_lab <- iCFCV(K=5,
                  treeNo=2000, 
                  iterationNo=100,
                  min.split.var=4, 
                  P_threshold=0.1, 
                  variable_type = "non-HD",
                  hdpct= 0.95,
                  HTE_P_cf.raw = 0.1
)
  
```
