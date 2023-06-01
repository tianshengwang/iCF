# iCF
iterative Causal Forest (iCF)

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
```

