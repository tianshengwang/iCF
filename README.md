# Iterative Causal Forest (iCF): A Novel Algorithm for Subgroup Identification
----------------------------------------------------------------------------------
The iCF algorithm identifies important subgroups with heterogeneous treatment effects without prior knowledge of treatment-covariate interactions

<img src = images/FIG1_19Mar2023.jpg width=1000>

**Citation**

Wang T, Keil AP, Kim S, Wyss R, Htoo PT, Funk MJ, Buse JB, Kosorok MR, Stürmer T. Iterative Causal Forest: A Novel Algorithm for Subgroup Identification. Am J Epidemiol. 2023 (In Press).

**R packages recommended**
```{r packages, include=FALSE}
library(MASS)
library(grf)
library(tidyverse)
library(rlang)
library(rlist)
library(plyr)
library(caret)
library(caTools)
library(listdtr)
library(randomForest)
library(ggplot2)
library(ggridges)
library(data.table)
library(grid)
library(broom)
library(rstatix)
library(DMwR)
library(knitr)
library(Rfast)
library(spaMM)
```
**installation**

There are three ways to install iCF.

First, download all the files in the 'R' folder and save them in your local directory (create a new folder named 'iCF'). Then, run the following R codes each time before executing iCF:
```{}
source("/local/iCF/best_tree_MSegar.R")
source("/local/iCF/iCF_TREE_build.R")
source("/local/iCF/iCF_PARENT_node.R")
source("/local/iCF/iCF_PRE_majority.R")
source("/local/iCF/iCF_MAJORITY_VOTE.R")
source("/local/iCF/iCF_SUBGROUP_DECISION.R")
source("/local/iCF/iCF_SG_PIPELINE.R")
source("/local/iCF/iCF_CV.R")
source("/local/iCF/iCF_SUBGROUP_ANALYSIS.R")
source("/local/iCF/iCF_iCF_grplasso.R")
source("/local/iCF/sim_Truth_tree.R")
source("/local/iCF/GG_toolbox.R")
```
Second, you can install it using devtools (this option is currently being worked on and will be available soon):
```{}
install.packages('devtools')
devtools::install_github("tianshengwang/iCF")
library(iCF)
```
Third, you can download the 'iCF_0.0.0.9000.tar.gz' file (this option is currently being worked on and will be available soon)..."

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
 W <- rbinom(nstudy,1,PS) 
#Two way interaction of W, X1, and X3
 Y = a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*W + rnorm(nstudy,0,1) + 0.4*W*X3 + 0.3*W*X1 + 0.4*W*X1*X3 + 0.2*X1*X3 
 dat <<- as.data.frame(cbind(W, Y, X1, X2, X3 ,X4, X5, X6, X7, X8, X9, X10)) 
``` 
**Run iCF**
***Run raw causal forest to predict outcome (Y.hat), propensity score (W.hat), and select variables***
```{}
 vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )
 X <- dat[,vars_forest]
 Y <- as.vector( as.numeric( dat[,"Y"] ) )
 W <- as.vector( as.numeric( dat[,"W"] ) )
 
 cf_raw_key.tr <- CF_RAW_key(dat, 1, "non-hd", hdpct=0.90) 
 Y.hat  <<- cf_raw_key.tr$Y.hat                 
 W.hat  <<- cf_raw_key.tr$W.hat                 
 HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    
 varimp_cf  <- cf_raw_key.tr$varimp_cf          
 selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx 

 GG_VI(varimp_cf, "Variable importance" )
 ```
 <img src = images/GG_VI_fig.png width=300>
 
 ***Tune leafsize to grow D2, D3, D4, and D5 causal forest***
 ```{}
#Specify the decimal position for continuous variables in the subgroup definition.
split_val_round_posi=0
#For real-world projects (not simulations where we know the truth), the truth is set as "Unknown".
truth.list <<- TRUTH("Unknown")
#Define categorical variables with more than two levels.
vars_catover2 <<- NA  
```
```{}
D2_MLS=MinLeafSizeTune(denominator=25, treeNo = 1000, iterationNo=100, "D2")
D2_MLS$depth_mean
D2_MLS$depth_gg
```
<img src = images/D2_MLS_tune.png width=400>

```{}
D3_MLS=MinLeafSizeTune(denominator=45, treeNo = 1000, iterationNo=100, "D3")
D3_MLS$depth_mean
D3_MLS$depth_gg
```
<img src = images/D3_MLS_tune.png width=350>

```{}
D4_MLS=MinLeafSizeTune(denominator=65, treeNo = 1000, iterationNo=100, "D4")
D4_MLS$depth_mean
D4_MLS$depth_gg
```
<img src = images/D4_MLS_tune.png width=350>

```{}
D5_MLS=MinLeafSizeTune(denominator=85, treeNo = 1000, iterationNo=100, "D5")
D5_MLS$depth_mean
D5_MLS$depth_gg
```
<img src = images/D5_MLS_tune.png width=350>

***Implement iCF***
```{}
leafsize <<- list(D5=85, D4=65, D3=45, D2=25)

iCFCV_B1000_i200 <- iCFCV(dat=dat,
                          K=5,
                          treeNo=1000, 
                          iterationNo=100,
                          min.split.var=4, 
                          P_threshold=0.1, 
                          variable_type = "non-HD",
                          hdpct= 0.95,
                          HTE_P_cf.raw = HTE_P_cf.raw) 

iCFCV_B1000_i200
```
If you have further questions or comments, please contact Dr. Tiansheng Wang: tianwang@unc.edu
