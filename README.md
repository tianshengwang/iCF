# Iterative Causal Forest (iCF): A Novel Algorithm for Subgroup Identification 

## iCF is free for nonprofit use

**Citation**

**Wang T, Keil AP, Kim S, Wyss R, Htoo PT, Funk MJ, Buse JB, Kosorok MR, Stürmer T. Iterative Causal Forest: A Novel Algorithm for Subgroup Identification. _Am J Epidemiol._ 2023 (In Press).**



The iCF algorithm (based on [causal forest](https://grf-labs.github.io/grf/articles/diagnostics.html)) identifies important subgroups with heterogeneous treatment effects without prior knowledge of treatment-covariate interactions

<img src = images/FIG1_19Mar2023.jpg width=1000>

**Iterative causal forest (iCF) algorithm**. CF, causal forest; RF, random forest; $\hat{Y}$, predicted outcome; $\hat{W}$, predicted propensity score; CV, cross-validation; Y*, predicted transformed outcome; $D_{tr}$, training data; $G_{D}$ subgroup decision; $g_{j}$, individual group of $G_{D}$. A. iCF workflow. If the homogeneity test from a raw casual forest is significant at P-value $\alpha$ = 0.1, then divide the data into 10 groups of equal size. The first group is treated as the testing set, the remaining 9 groups as the training set. Repeat the procedure 10 times. B. Step 3 of iCF workflow (obtain a family of subgroup decisions from a pruned iCF on the training set and a family of transformed outcome models based on subgroup decisions $G_{D}$ on the test stet). With a cross-validation approach, each $G_{D}$ selected will be noted as $G_{D_f}$, where f denotes the fold of cross-validation. If $G_{D_f}$ varies across training sets, then perform plurality vote across f subgroup decisions to select the most popular one and use its corresponding model. Then, average the error across validation sets for each model from $G_{2}$, $G_{3}$, $G_{4}$, and $G_{5}$, individually, and the $G_{D}$ with a corresponding model with the smallest cross-validated error is selected as the final subgroup decision $G_{iCF}$



**1. R packages recommended**
```{r packages, include=FALSE}
library(MASS)
library(grf)
library(tidyverse)
library(rlang)
library(rlist)
library(plyr)
library(caret)
library(caTools)
library(randomForest)
library(data.table)
library(grid)
library(broom)
library(rstatix)
library(knitr)
library(Rfast)
library(ggplot2)
library(ggridges)
library(glmnet)
```
**2. Installation**

There are three ways to install iCF.

First, download all the R programs in the [R folder](https://github.com/tianshengwang/iCF/tree/iCF/R) and save them in your local directory (create a new folder named 'iCF'). Then, run the following R codes each time before executing iCF:
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
source("/local/iCF/iCF_grplasso.R")
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

```{}
install.packages("iCF_0.0.0.9000.tar.gz", repos = NULL, type ="source")
library(iCF)
```

**3. Data Simulation**

We used a simulation setup similar to the one initially described by [Setoguchi et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2905676/) and modified it by incorporating treatment-covariate interactions to model heterogeneous treatment effects.

<img src = images/SimulationData.png width=600>

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

#4 standard normal: X2, X4, X7, X10;
  X2 <- rnorm(nstudy, mean=0, sd=1)  
  X4 <- rnorm(nstudy, mean=0, sd=1)
  X7 <- rnorm(nstudy, mean=0, sd=1)    
  X10 <- rnorm(nstudy, mean=0, sd=1)

#6 binary variables: X1, X3, X5, X6, X8, X9;
  X1 <- rnorm(nstudy, mean=0, sd=1)
  X1 <- ifelse(X1 > mean(X1), 1, 0)

  X3 <- rnorm(nstudy, mean=0, sd=1)
  X3 <- ifelse(X3 > mean(X3), 1, 0)

  X5 <- F.sample.cor(X1, 0.2)
  X5 <- ifelse(X5 > mean(X5), 1, 0)

  X6 <- F.sample.cor(X2, 0.9)
  X6 <- ifelse(X6 > mean(X6), 1, 0)

  X8 <- F.sample.cor(X3, 0.2)
  X8 <- ifelse(X8 > mean(X8), 1, 0)

  X9 <- F.sample.cor(X4, 0.9)
  X9 <- ifelse(X9 > mean(X9), 1, 0)

#4 confounders: X1, X2, X3, X4;
#3 outcome predictors: X8, X9, X10;
#1 instrumental X7;
 
 PS <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7) ))^-1 #true propensity score
 W <- rbinom(nstudy,1,PS) 
#3-way interaction of W, X1, and X3 (W:X1:X3)
 Y = a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*W + rnorm(nstudy,0,1) + 0.4*W*X3 + 0.3*W*X1 + 0.4*W*X1*X3 + 0.2*X1*X3
```
Make a dataset in this format to be run by iCF: the 1st column is treatment **W**, 2nd column is outcome **Y**, and remaining columns features **X**.
```{}
 Train <<- as.data.frame(cbind(W, Y, X1, X2, X3 ,X4, X5, X6, X7, X8, X9, X10)) 
``` 
**4. Run iCF on simulated data**

***Step 1. Run raw causal forest to predict outcome (Y.hat), propensity score (W.hat), and select variables***
```{}
 vars_forest = colnames( Train %>% dplyr::select(-c("Y", "W" ))  )
 X <- Train[,vars_forest]
 Y <- as.vector( as.numeric( Train[,"Y"] ) )
 W <- as.vector( as.numeric( Train[,"W"] ) )
 cf_raw_key.tr <- CF_RAW_key(Train, 1, "non-hd", hdpct=0.90) 
 Y.hat  <<- cf_raw_key.tr$Y.hat                 
 W.hat  <<- cf_raw_key.tr$W.hat                 
 HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    
 varimp_cf  <- cf_raw_key.tr$varimp_cf          
 selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx
 X[,selected_cf.idx]%>%colnames()
 GG_VI(varimp_cf, "Variable importance", colnames(X))
 ```
 <img src = images/GG_VI_fig.png width=400>
 
 ***Step 2: Tune the minimum leaf size (MLS) for D2, D3, D4, and D5 to ensure that the majority of the best trees from causal forests grown with these MLS have depths of 2, 3, 4, and 5, respectively.***
 ```{}
#Specify the decimal position for continuous variables in the subgroup definition.
split_val_round_posi=0
#Define categorical variables with more than two levels. If no such variables, let it be NA:
vars_catover2 <<- NA  
```
```{}
D2_MLS=MinLeafSizeTune(Train, denominator=25, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D2", "steelblue1")
D2_MLS$depth_mean
D2_MLS$depth_gg
```
<img src = images/D2_MLS_tune.png width=400>

Notably, if you got this message "_Error: Can't subset columns that don't exist. x Column `parent_sign` doesn't exist._", it suggests the denominator used for developing is too small, leading to a too large MLS for D2 forest so that the tree does not even split (the node does not have a parent node). In this scenario, increasing the denominator will solve the problem.


```{}
D3_MLS=MinLeafSizeTune(Train, denominator=45, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D3", "steelblue1")
D3_MLS$depth_mean
D3_MLS$depth_gg
```
<img src = images/D3_MLS_tune.png width=400>

```{}
D4_MLS=MinLeafSizeTune(Train,denominator=65, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D4", "steelblue1")
D4_MLS$depth_mean
D4_MLS$depth_gg
```
<img src = images/D4_MLS_tune.png width=400>

```{}
D5_MLS=MinLeafSizeTune(Train, denominator=85, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D5", "steelblue1")
D5_MLS$depth_mean
D5_MLS$depth_gg
```
<img src = images/D5_MLS_tune.png width=400>

*Note that despite using a smaller MLS, the best trees from causal forests do not grow deeper due to the presence of a strong three-way interaction (W:X1:X3) in the simulated data set. In such scenarios, we can proceed with implementing iCF when tuning MLS does not affect the depth of the best tree.* 

***Step 3. Implement iCF on simulated dataset***
```{}
leafsize <<- list(D5=D5_MLS$denominator, D4=D4_MLS$denominator, D3=D3_MLS$denominator, D2=D2_MLS$denominator)

iCFCV_B1000_i200_sim <- iCFCV(dat=Train, K=5, treeNo=200, iterationNo=10, min.split.var=4,
                              split_val_round_posi=0, P_threshold=0.5, variable_type = "non-HD", 
                              hdpct= 0.95, HTE_P_cf.raw = HTE_P_cf.raw) 

iCFCV_B1000_i200_sim
```


**5. Run iCF on real-world data**

We compared the **two-year risk difference** of hospitalized heart failure (HHF) of initiating any sodium-glucose cotransporter-2 inhibitors (SGLT2i) versus glucagon-like peptide-1 receptor agonists (GLP1RA) using a 20% random sample of all fee-for-service U.S. Medicare beneficiaries who had parts A (inpatient), B (outpatient physician services), and D (dispensed prescription drugs) coverage for at least one month from January 2012 to December 2017. The details of the cohort were described previously by [Htoo et al.](https://www.ahajournals.org/doi/full/10.1161/JAHA.121.022376) and are available in the mehtod paper (Wang et al.) 

Note the outcome Y should be on **risk difference** scale. For example, only 2 (Patient 1 & 4) out of the following 5 patients had heart failure hospitalization (black dot) events in a 2-year period.


<img src = images/RD_scale_outcome.png width=800>


Make a dataset in this format to be run by iCF: the 1st column is treatment **W**, 2nd column is outcome **Y**, and remaining columns features **X**. And the 1st and 2nd columns should be named following this pattern, i.e., 'W' for treatment and 'Y' for outcome.

```{}
load("HHF_SGLTvGLP_iCFCV.RData")

Train <-  hfp_2yr_all_sgltvglp %>% 
          dplyr::filter(FillDate2 <= 21184 -365*2 &  #31DEC2017
                        IndexDate >= 19449           #1APR2013
                        ) %>%
          dplyr::mutate(Y = hfp_2yr_2yr,
                        age = as.numeric(cut(age, c(65,70,75,80,85,Inf) ,
                                             labels=c("65<age<=70 ","70<age<=75","75<age<=80","80<age<=85", "age>85")))) %>% 
          dplyr::rename(W=SGLT) %>%
          dplyr::select(Y, W, age, race2, sex,baselinecvd, baselinechf, 
                        dplyr::starts_with("bl_")) %>% 
                      as.data.frame.matrix() %>%
                      mutate(sex=as.numeric(sex))
```
***Step 1. Run raw causal forest to predict outcome (Y.hat), propensity score (W.hat), and select variables***

```{}
#Define categorical variables with more than two levels.
vars_catover2 <- c("race2", "age", "bl_HOSP", "bl_HOSPDAYS", "bl_ED", "bl_ERDM", "bl_outpt", "bl_outptdm")

vars_forest = colnames( Train %>% dplyr::select(-c("Y", "W"))  ) 
X <<- Train[,vars_forest]
Y <<- as.vector( as.numeric( Train[,"Y"] ) )
W <<- as.vector( as.numeric( Train[,"W"] ) )

cf_raw_key.tr <- CF_RAW_key(Train, min.split.var=4, variable_type="non-hd", hdpct=0.95)    
Y.hat  <<- cf_raw_key.tr$Y.hat
W.hat  <<- cf_raw_key.tr$W.hat
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw
varimp_cf  <- cf_raw_key.tr$varimp_cf
selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx
time_rawCF <- cf_raw_key.tr$time_rawCF
X[,selected_cf.idx]
GG_VI(varimp_cf, 'Variable Importance for SGLT2i vs GLP1RA cohort for HFF', colnames(X))
 ```
 <img src = images/VI_HHF2y_nolabel.png width=800>
 
 ***Step 2: Tune the MLS for D2, D3, D4, and D5 to ensure that the majority of the best trees from causal forests grown with these MLS have depths of 2, 3, 4, and 5, respectively.***

```{}
D2_MLS=MinLeafSizeTune(dat=Train, denominator=25, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D2", "#62C6F2")
D2_MLS$depth_mean
D2_MLS$depth_gg
```
<img src = images/D2_MLS_tune_rwd.png width=400>

```{}
D3_MLS=MinLeafSizeTune(dat=Train, denominator=45, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D3", "#62C6F2")
D3_MLS$depth_mean
D3_MLS$depth_gg
```
<img src = images/D3_MLS_tune_rwd.png width=400>

```{}
D4_MLS=MinLeafSizeTune(dat=Train, denominator=65, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D4", "#62C6F2")
D4_MLS$depth_mean
D4_MLS$depth_gg
```
<img src = images/D4_MLS_tune_rwd.png width=400>

```{}
D5_MLS=MinLeafSizeTune(dat=Train, denominator=85, treeNo = 200, iterationNo=10, split_val_round_posi=0, "D5", "#62C6F2")
D5_MLS$depth_mean
D5_MLS$depth_gg
```
<img src = images/D5_MLS_tune_rwd.png width=400>

***Step 3. Implement iCF on Medicare SGLT2i vs GLP1RA new user cohort***
```{}
leafsize <<- list(D5=D5_MLS$denominator, D4=D4_MLS$denominator, D3=D3_MLS$denominator, D2=D2_MLS$denominator)

iCFCV_B1000_i200_rwd <- iCFCV(dat=Train,K=5, treeNo=200, iterationNo=25, min.split.var=4,
                              split_val_round_posi=0, P_threshold=0.5, variable_type = "non-HD", 
                              hdpct= 0.95, HTE_P_cf.raw = HTE_P_cf.raw) 

iCFCV_B1000_i200_rwd
```

If you have further questions or comments, please contact Dr. Tiansheng Wang: tianwang@unc.edu
