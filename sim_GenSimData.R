#############################################################################################
# DATA PREPARE 
# Author: Tiansheng Wang  (Simulation data part based on Brian K Lee's simulation's R code )
# Last update date: 9/6/2020
# Version: 0.1         
#############################################################################################

#' This function returns the inverse logit of a given numeric value, i.e., translate it into probabilities between 0 and 1. 
#' @param m the model of event Y (binary outcome), m = alpha + beta*X
#' 
#' @return The probability of Y
#'
#' @export
#'
inv.logit <- function(m){
  p=exp(m)/(1+exp(m)) 
  return(p)
}




#' This function generate simulated data, may adjust parameter to check different effect size
#' @param scenario scenario for propensity score model
#' @param intTRUE scenario for HTE
#' @param a8 coefficient for W:X1
#' @param a9 coefficient for W:X3
#' @param a10 coefficient for X1:X3
#' @param a11 coefficient for W:X1:X3 
#' @param a12 coefficient for W:X8
#' @param a14 coefficient for X3:X8
#' @param a15 coefficient for W:X1:X8
#' @param a16 coefficient for W:X3:X8
#' @param a17 coefficient for X1:X3:X8
#' @param a18 coefficient for W:X1:X3:X8
#' @param a21 coefficient for W:X2
#' @param a22 coefficient for X2:X3
#' @param a23 coefficient for W:X2:X3
#' @param a2s2 coefficient for I(X2^2) 
#' @param a4s2 coefficient for I(X4^2)
#' @param a1i3 coefficient for non-additivity X1:X3, if a1i3!=0 then a10=0
#' @param a2i4 coefficient for non-additivity X2:X4 
#' @param a2i3 coefficient for non-additivity X2:X3, if a2i3!=0 then a22=0
#' @param a3i4 coefficient for non-additivity X3:X4
#' @param nstudy sample size

#' 
#' @return The subsetted list that all are identical to reference list.
#'
#' @export
#
GEN_SIM_DATA_HTE <- function(scenario, intTRUE, outcome_type, a8,	  a9,	  a10,  a11,   a12, a12_, a13, a13_, a14, a14_,  a15, a15_, a16, a16_, a17, a17_, a18, a18_, a21, a22, a23, a2s2, a4s2, a1i3, a2i4, a2i3, a3i4, nstudy, seed=NULL) {
 
   set.seed(seed)
  
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
    #cat("Sample corr = ", cor(x, y), "\n")
    return(Y)
  }
  
#  set.seed()
  #scenario<- scenar		#scenario<-c(rep(scenar,nstudy))
  #----------------------------------------
  #generating covariates (columns)
  #----------------------------------------
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
  #~~ dichotomize variables (will attenuate correlations above)
  X1 <- ifelse(X1 > mean(X1), 1, 0)
  X3 <- ifelse(X3 > mean(X3), 1, 0)
  X5 <- ifelse(X5 > mean(X5), 1, 0)
  X6 <- ifelse(X6 > mean(X6), 1, 0)
  X8 <- ifelse(X8 > mean(X8), 1, 0)
  X9 <- ifelse(X9 > mean(X9), 1, 0)
  #~~ scenarios for data generation models
  # A: model with additivity and linearity
  # B: mild non-linearity
  # C: moderate non-linearity
  # D: mild non-additivity
  # E: mild non-additivity and non-linearity
  # F: moderate non-additivity
  # G: moderate non-additivity and non-linearity
  # binary exposure modeling
  
  if (scenario == "A") {
    z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7) ))^-1
    #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7)			#probit model
    testing<- 1
  } else
    if (scenario == "B") {
      z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7 + b2*X2*X2) ) )^-1
      #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7 + b2*w2*w2)			#probit model
      testing<-2
    } else
      if (scenario == "C") {
        z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7 + b2*X2*X2 +b4*X4*X4 + b7*X7*X7) ) )^-1
        #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7 + b2*w2*w2 +b4*w4*w4 + b7*w7*w7)			#probit model
        testing<- 3
      } else
        if (scenario == "D") {
          z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7 + b1*0.5*X1*X3 + b2*0.7*X2*X4 + b4*0.5*X4*X5 + b5*0.5*X5*X6) ) )^-1
          #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7 + b1*0.5*w1*w3 + b2*0.7*w2*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6)			#probit model
          testing<- 4
        } else
          if (scenario == "E") {
            z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7
                                      + b2*X2*X2 + b1*0.5*X1*X3 + b2*0.7*X2*X4 + b4*0.5*X4*X5 + b5*0.5*X5*X6) ) )^-1
            #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
            #			+ b2*w2*w2 + b1*0.5*w1*w3 + b2*0.7*w2*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6)			#probit model
            testing<- 5
          } else
            if (scenario == "F") {
              z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7
                                        + b1*0.5*X1*X3 + b2*0.7*X2*X4 + b3*0.5*X3*X5 + b4*0.7*X4*X6 + b5*0.5*X5*X7
                                        + b1*0.5*X1*X6 + b2*0.7*X2*X3 + b3*0.5*X3*X4 + b4*0.5*X4*X5 + b5*0.5*X5*X6) ) )^-1
              #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
              #			+ b1*0.5*w1*w3 + b2*0.7*w2*w4 + b3*0.5*w3*w5 + b4*0.7*w4*w6 + b5*0.5*w5*w7
              #			+ b1*0.5*w1*w6 + b2*0.7*w2*w3 + b3*0.5*w3*w4 + b4*0.5*w4*w5 + b5*0.5*w5*w6)			#probit model
              testing<- 6
            } else
            { #scenario G
              z.a_trueps <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7
                                        + b2*X2*X2 + b4*X4*X4 + b7*X7*X7 + b1*0.5*X1*X3 + b2*0.7*X2*X4 +b3*0.5*X3*X5
                                        + b4*0.7*X4*X6 + b5*0.5*X5*X7 + b1*0.5*X1*X6 + b2*0.7*X2*X3 + b3*0.5*X3*X4
                                        + b4*0.5*X4*X5 + b5*0.5*X5*X6) ) )^-1
              #z.a_trueps <- pnorm(b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7
              #			+ b2*w2*w2 + b4*w4*w4 + b7*w7*w7 + b1*0.5*w1*w3 + b2*0.7*w2*w4 +b3*0.5*w3*w5
              #			+ b4*0.7*w4*w6 + b5*0.5*w5*w7 + b1*0.5*w1*w6 + b2*0.7*w2*w3 + b3*0.5*w3*w4
              #			+ b4*0.5*w4*w5 + b5*0.5*w5*w6)									#probit model following the example of Brookhart et al.
              testing<- 7
            }
  
  # probability of exposure: random number betw 0 and 1
  # if estimated true ps > prob.exposure, than received exposure (z.a=1)
  #prob.exposure <- runif(nstudy)
  #z.a <- ifelse(z.a_trueps > prob.exposure, 1, 0)
  a<- rbinom(nstudy,1,z.a_trueps)
  
 # X4w=X8

  X2_con = X2
  X2df <- as.data.frame(cbind(X2, X2_con)) %>% 
    mutate(X2_l3 = ifelse(X2 > quantile(X2_con,0.75) , 2,
                          ifelse(X2 <=  quantile(X2_con,0.75) & X2 > quantile(X2_con,0.25), 1,  0 )),
           
           X2_l4 = ifelse(X2 > quantile(X2_con,0.75) , 3,
                          ifelse(X2 <=  quantile(X2_con,0.75) & X2 > quantile(X2_con,0.5), 2,
                                 ifelse(X2 <=  quantile(X2_con,0.5) & X2 > quantile(X2_con, 0.25), 1, 0)))
    )
  if ( stringr::str_sub(intTRUE, -2, -1) == "o4"){
    X2 <- X2df$X2_l4  
    
  } else if ( stringr::str_sub(intTRUE, -2, -1) == "o3"){
    X2 <- X2df$X2_l3
    
  }
  #_______
#___________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________	  
  # continuous outcome modeling, be careful!!!don't separate by 3 lines, express y in one line!!!
  Y_main <- a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*a + rnorm(nstudy,0,1)
  
  Y_binary_int <- a9*a*X3 + a8*a*X1 + a10*X1*X3  + a11*a*X1*X3 + a12*a*X8 + a13*X1*X8 + a14*X3*X8  + a15*a*X1*X8  + a16*a*X3*X8  + a17*X1*X3*X8 + a18*a*X1*X3*X8 
  
  #----------------------------------------------------------------
  #X2 has 4 levels: 0, 1, 2, 3 or X2 is a continuous variable
  #----------------------------------------------------------------
  if ( stringr::str_sub(intTRUE, -2, -1) == "o4"|stringr::str_sub(intTRUE, -2, -1) == "o3"){
    
  Y_nonbinary_int <- a21*a*(X2>1)   + a22*(X2>1)*X3   + a23*a*(X2>1)*X3   + a12_*a*(X2>1)   + a13_*X1*(X2>1)   + a14_*X3*(X2>1)   + a15_*a*X1*(X2>1)   + a16_*a*X3*(X2>1)   + a17_*X1*X3*(X2>1)   + a18_*a*X1*X3*(X2>1)
  
  } else if ( stringr::str_sub(intTRUE, -1, -1) == "l"){
  
    #update Truth_tree.R program accordingly to be consistent!
  Y_nonbinary_int <- a21*a*(X2>0.5) + a22*(X2>0.5)*X3 + a23*a*(X2>0.5)*X3 + a12_*a*(X2>0.5) + a13_*X1*(X2>0.5) + a14_*X3*(X2>0.5) + a15_*a*X1*(X2>0.5) + a16_*a*X3*(X2>0.5) + a17_*X1*X3*(X2>0.5) + a18_*a*X1*X3*(X2>0.5)
  
  }  else {
    
  Y_nonbinary_int=0
  }
    
  
  #non-linearity: add suare term
  Y_nonline <- a2s2*X2*X2 + a4s2*X4*X4
  #non-additivity: add interaction term
  Y_nonadit <- a1*a1i3*X1*X3 + a2*a2i4*X2*X4 + a2*a2i3*X2*X3 + a3*a3i4*X3*X4
  #sum to get final outcome model
 Y_hte =  Y_binary_int + Y_nonbinary_int + Y_nonline + Y_nonadit
 
         if (outcome_type =="continuous") {
  #continuous outcome
  Y1          =             Y_main +  Y_hte
  dat <<- as.data.frame(cbind(X1, X2, X3 ,X4, X5, X6, X7, X8, X9, X10, z.a_trueps, a, Y1))
  
  } else if (outcome_type =="binary")     {
  #binary outcome
#interacept a0=0.1 to replace  original value - 3.85, if negative, e.g. even -0.1, it is very likely to get a negative probability
#risk difference scale: "everything with causal forests works the same for both continuous and binary outcomes. If you have binary outcomes, then the CATE is estimated on the "difference in probabilities" scale". https://github.com/grf-labs/grf/issues/238
#thus logistic sigmoid function (1 + exp(-x))^-1 only applied to covariates X, leave both treatment and interaction terms out
# (1 + exp(-0.35*x))^-1: add 0.35 to make sigmoid.x in an appropriate range so that prob.y 
  #sigmoid.x <-   (1+ exp( 0.1*(  a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 ) ) )^-1
  sigmoid.x <-   (1+ exp( -0.1*(  a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 ) ) )^-1
  
  hist(sigmoid.x); range(sigmoid.x)
  #hist(Y_hte);     range(Y_hte)
  #hist(g1*a );     range(g1*a )
  #hist( g1*a +  sigmoid.x +  Y_hte);     range( g1*a +  sigmoid.x +  Y_hte)
  #hist( 0.01 + g1*a + sigmoid.x + Y_hte);     range(0.01 + g1*a + sigmoid.x + Y_hte)
  #hist( 0.01 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*a  + Y_hte);     range(0.01 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*a + Y_hte)
  
#g1=-0.4;  new g1=-0.12
#a9=0.4;   new a9=0.12
  prob.y <- 0.01 +  g1*a +  sigmoid.x +  Y_hte +  rnorm(nstudy,0,0.02)
  hist(prob.y);   range(prob.y)
  mean(prob.y)
  
  Y1 = rbinom (nstudy, 1, prob.y)
  mean(Y1); table(Y1)
  dat <<- as.data.frame(cbind(X1, X2, X3 ,X4, X5, X6, X7, X8, X9, X10, z.a_trueps, a, Y1, prob.y))
  
  }
  #notice: "<<-" is what saves "dat" in your environment

  return(dat)
  
}

########################################

