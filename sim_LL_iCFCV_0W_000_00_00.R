####################################################################
#       HTE () w/o UC
#  Author: Tiansheng Wang, 
# Last update date: 9/23/2020
# Version: 0.1   
####################################################################
rm(list = ls())
library(MASS)
library("grf")
library(tidyverse)
library(rlang)
library(rlist)
library(plyr)
library(policytree)
library(caret)
library(caTools)
library(listdtr)
library(randomForest)
library("rlist")
library(ggplot2)
library(ggridges)
theme_set(theme_ridges())
library(DiagrammeRsvg)
library("DiagrammeR")
library(data.table)
library(grid)
library(cobalt)
library(broom)
library(optmatch)
#library(stima)
library(personalized)
library(FindIt)
library(grplasso)
library(rstatix)
library(ggpubr)
library(DMwR)
library(aVirtualTwins)
library(knitr)
library(Rfast)
library("e1071")
library("rcdd")
library("spaMM")
library(DynTxRegime)

setwd("/nas/longleaf/home/tianwang/iCF/output/")

#keep this order of source
source("/nas/longleaf/home/tianwang/iCF/programs/macros/GG_toolbox.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/sim_GenSimData.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/best_tree_MSegar.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_TREE_build.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_PARENT_node.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_PRE_majority.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/sim_Truth_tree.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_MAJORITY_VOTE.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_SUBGROUP_DECISION.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_SG_PIPELINE.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_SUBGROUP_ANALYSIS.R")
#Interaction trees
source("/nas/longleaf/home/tianwang/iCF/programs/macros/Functions-IT.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/Functions-RFIT.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/DTR_personalized.R")
source("/nas/longleaf/home/tianwang/iCF/programs/macros/iCF_grplasso.R")


###############################################################################################
# causal forest wrapper function for easy-input of different parameter
###############################################################################################
#doesn't work if run from source for real simulation, which shows error "Validate_X, can't find object X", thus run here!
CF <- function(depth, treeNo, tunepara){
  CausalForest <-grf::causal_forest(X[, selected_cf.idx],     
                                    Y, 
                                    W,  
                                    Y.hat ,
                                    W.hat ,
                                    num.trees = treeNo,
                                    sample.weights = NULL,
                                    clusters = NULL,
                                    equalize.cluster.weights = FALSE,
                                    sample.fraction = 0.5,
                                    mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                                    min.node.size  = round(nrow(Train)/depth, 0), 
                                    honesty = TRUE,
                                    honesty.fraction = 0.5,
                                    honesty.prune.leaves = TRUE,
                                    alpha = 0.05,
                                    imbalance.penalty = 0,
                                    stabilize.splits = TRUE,
                                    ci.group.size  = 2,
                                    tune.parameters = tunepara, #only "none" can produce simple forest,
                                    tune.num.trees = 200,
                                    tune.num.reps = 50,
                                    tune.num.draws = 1000,
                                    compute.oob.predictions = TRUE,
                                    #orthog.boosting = FALSE,#Deprecated and unused after version 1.0.4.
                                    num.threads = NULL,
                                    seed =  runif(1, 0, .Machine$integer.max))
  return(CausalForest)
}

###############################################################################################
#                      Iterative multidepth CF  + Majority Vote                               #
###############################################################################################
iCF <- function( depth, treeNo, iterationNo, Ntrain, tree_depth, split_val_round_posi){
  besttreelist = list() #making an empty list
  besttreelist_L = list () #for list format (rather than df format) of best trees
  splitfreqlist = list()
  treeBlist <- list()
  #treeBlist_L <- list()
  mvtreelist = list()
  # mvtreelist_L = list()
  for (k in 1:iterationNo) { 
    cf <- CF(depth, # =as.numeric(sample(leafsize, 1)), 
             treeNo, "none")
    #cf <- CF(depth, treeNo, "none")
    #-----------------------------------------------------------------
    #get P-value of omnibus test for HTE presence 
    HTE_P_cf <- grf::test_calibration(cf)[8]
    #-----------------------------------------------------------------
    #Prepare for Split Frequency
    max.depth_para= as.numeric(stringr::str_sub(tree_depth,-1,-1))
    #split frequency for whole forest, not frequency of the "depth"!
    freqs <- grf::split_frequencies(cf, max.depth = max.depth_para)
    d <- data.frame(freqs)
    real_index <- colnames(X[, selected_cf.idx])
    #rename with real index 
    data.table::setnames(d, old = c(colnames(d)), new = c(real_index))
    d$k <- k
    splitfreqlist[[k]] <- d # add it to your list
    #-----------------------------------------------------------------
    #############################-----------------------------------------------------------------
    #BEST TREE SELECTION
    #############################-----------------------------------------------------------------
    best_tree_info<-find_best_tree(cf, "causal")
    best_tree_info$best_tree
    # Plot trees
    tree.plot = plot(grf::get_tree(cf, best_tree_info$best_tree))
    tree.plot
    best.tree <- grf::get_tree(cf, best_tree_info$best_tree)
    best.tree
    besttreelist_L[[k]] <- best.tree # add it to your list
    #-----------------------------------------------------------------
    besttree <-   GET_TREE_DF(length(best.tree$nodes), best.tree, split_val_round_posi)
    besttree$k <- k  # maybe you want to keep track of which iteration produced it?
    besttree$HTE_P_cf <- HTE_P_cf  # record the omnibus HTE test P value
    besttreelist[[k]] <- besttree # add it to your list
  }
  
  all_list <- rlist::list.zip(besttreelist, besttreelist_L, splitfreqlist)
  return(all_list) 
}


###############################################################################################
###############################################################################################
#         ONLY ONE multidepth CF  + Majority Vote,i.e. without best tree !!!                  #
###############################################################################################
###############################################################################################

oneCF <- function(depth, treeNo, tunepara, tree_depth, split_val_round_posi){
  treeBlist  = list()
  besttreelist = list() #making an empty list
  besttreelist_L = list () #for list format (rather than df format) of best trees
  cf <- CF(depth, treeNo, tunepara)
  #-----------------------------------------------------------------
  #get P-value of omnibus test for HTE presence 
  HTE_P_cf <- grf::test_calibration(cf)[8]
  #-----------------------------------------------------------------
  #Prepare for Split Frequency
  max.depth_para= as.numeric(stringr::str_sub(tree_depth,-1,-1))
  #split frequency for whole forest, not frequency of the "depth"!
  freqs <- grf::split_frequencies(cf, max.depth = max.depth_para)
  d <- data.frame(freqs)
  real_index <- colnames(X[, selected_cf.idx])
  #rename with real index 
  data.table::setnames(d, old = c(colnames(d)), new = c(real_index))
  #-----------------------------------------------------------------
  for (b in 1:treeNo) { 
    tree_b_l <- grf::get_tree(cf, b)
    tree_b <-   GET_TREE_DF(length(tree_b_l$nodes), tree_b_l, split_val_round_posi)
    tree_b$b <- b 
    tree_b$HTE_P_cf <- HTE_P_cf
    treeBlist[[b]] <- tree_b
  }
  #############################-----------------------------------------------------------------
  #BEST TREE SELECTION
  #############################-----------------------------------------------------------------
  best_tree_info<-find_best_tree(cf, "causal")
  best_tree_info$best_tree
  # Plot trees
  tree.plot = plot(grf::get_tree(cf, best_tree_info$best_tree))
  tree.plot
  best.tree <- grf::get_tree(cf, best_tree_info$best_tree)
  best.tree
  besttreelist_L[[1]] <- best.tree # add it to your list
  #-----------------------------------------------------------------
  besttree <-   GET_TREE_DF(length(best.tree$nodes), best.tree, split_val_round_posi)
  besttree$k <- 1  # maybe you want to keep track of which iteration produced it?
  besttree$HTE_P_cf <- HTE_P_cf  # record the omnibus HTE test P value
  besttreelist[[1]] <- besttree # add it to your list
  #############################-----------------------------------------------------------------
  return(list(besttreelist=besttreelist, besttreelist_L=besttreelist_L, treeBlist=treeBlist, cf=cf, HTE_P_cf=HTE_P_cf, d=d))
}
#source("/nas/longleaf/home/tianwang/iCF/programs/MAIN_miCF.R")

# function: generate continuous random variable correlated to variable x by rho
# invoked by the "F.generate" function
# Parameters -
# x - data vector
# rho - correlation coefficient
# Returns -
# a correlated data vector of the same length as x
# split_val_round_posi=1, if =0, i.e., rounding to integer, will lead to a subgroup with 0 observation, e.g. X7>0 & X7<=0 & X10<=0   
#---------------------


args = commandArgs(trailingOnly = TRUE)
rep = as.numeric(args[1])
set.seed(rep)






sim.function<- function(scenario, intTRUE, UC_indi, outcome_type, a8,	  a9,	  a10,  a11,  a12, a12_, a13, a13_, a14, a14_,  a15, a15_, a16, a16_, a17, a17_, a18, a18_,  a21, a22, a23, a2s2, a4s2, a1i3, a2i4, a2i3, a3i4, treeNo, iterationNo, nstart, nsim, nstudy, split_val_round_posi){
  
	for(j in nstart:nsim){

	  
	  
	  
	  intTRUE <<- intTRUE	  
split_val_round_posi <<- split_val_round_posi
	  #------------------------------------------------------------------		
	  #Get truth
	  truth.list <<- TRUTH(intTRUE)
	  #truth 
	  tree_true <<- truth.list$tree_true
	  #truth_r
	  tree_true_r <<- truth.list$tree_true_r
	  #truth_subgroup
	  tree_true_subgroup <<- truth.list$tree_true_subgroup
	  #truth_N1, 
	  tree_true_N1 <<- truth.list$tree_true_N1
	  #truth_N1_r, 
	  tree_true_N1_r <<- truth.list$tree_true_N1_r
	  #truth_N123
	  tree_true_N123 <<- truth.list$tree_true_N123
	  #truth_N123_r
	  tree_true_N123_r <<- truth.list$tree_true_N123_r
	  #truth description
	  truth_description <<- truth.list$truth_description
	  #truth INT
	  truth_INT <<- truth.list$truth_INT
	  #muticategorical variables, replace with vector of variable names for Medicare data!!!
	  vars_catover2 <<- truth.list$vars_catover2
	  #true subgroup decision to interaction (simplied, unified INT expression to assess performance)
	  true_SG2INT <<- SG2INT(tree_true_subgroup)
	  
	 
	  #------------------------------------------------------------------	
	  ######################################----------------------------------------------
	  ######################################----------------------------------------------
	  # I. Prepare for simulation data
	  ######################################----------------------------------------------
	  ######################################---------------------------------------------
	 dat00 <- GEN_SIM_DATA_HTE(scenario, intTRUE, outcome_type, a8,	  a9,	  a10,  a11,   a12, a12_, a13, a13_, a14, a14_,  a15, a15_, a16, a16_, a17, a17_, a18, a18_, a21, a22, a23, a2s2, a4s2, a1i3, a2i4, a2i3, a3i4, nstudy)
	
	 
	  #---------------------------------------------------------------------------------------------
	  ###       SECTION 2 RANDOMLY DIVIDED INTO 1/2 FOR TRAINING, 1/2 FOR TESTING
	  #---------------------------------------------------------------------------------------------
	  #https://www.listendata.com/2015/02/splitting-data-into-training-and-test.html
	  #for continuous outcome, no need split by prevalance of "binary outcome"
	  #trainIndex <- createDataPartition(sim$Y, p = .5, list = FALSE,times = 1)
	  #flip the sign to make Y>0 for 2 reasons:
	  #1) some ML based PM methods have default setting that the larger Y the better outcome
	  
	  if (outcome_type == "continuous"){
	  dat00 <- dat00 %>% 
	           dplyr:: mutate(Y= Y1 ) %>% #deactivate for real-world data 
	           dplyr:: rename(W=a) %>%  #deactivate for real-world data 
	           dplyr::select(-c(z.a_trueps, Y1)) %>%  #deactivate for real-world data 
	           dplyr::select(W, everything()) %>% 
	           dplyr::select(Y, everything())
	  } else if (outcome_type == "binary"){
	   
	    if(intTRUE !="Unknown"){
	     # outcome & %
	    table(dat00$Y1); prop.table(table(dat00$Y1))
	    # treatment & %
	    table(dat00$a); prop.table(table(dat00$a))
	    # outcome by treatment and %
	    table(dat00$Y1, dat00$a); prop.table(table(dat00$Y1, dat00$a))
	    prob.y <- dat00$prob.y
	    mean(prob.y); hist(prob.y); range(prob.y)
	    
	  dat00 <- dat00 %>% dplyr::select(-prob.y) %>%
	          dplyr:: mutate(Y= ifelse(Y1==0, 0, 1) ) %>% 
	          dplyr:: rename(W=a) %>% 
	          dplyr::select(-c(z.a_trueps, Y1)) %>% 
	          dplyr::select(W, everything()) %>% 
	          dplyr::select(Y, everything())
	    } else {
	  dat00 <- dat00 %>% dplyr::select(W, everything()) %>% 
	                    dplyr::select(Y, everything())  
	    }
	  
	  }

	  ############################################################
	  # Emperical study: where to load Medicar data "dat0" !!!
	  ############################################################

	  set.seed(20160413)
	  trainIndex = caret::createDataPartition(y=dat00[,1], p=0.8, list = FALSE)
	  
	  Train <<- dat00[ trainIndex,]
	  Test <<- dat00[-trainIndex,]
	  Ntrain <<- nrow(Train)
	  
	  ID <-1:nrow(dat00)
	  testIndex <- setdiff(ID, as.vector(trainIndex))
	  Train_ID <<- cbind(Train, as.vector(trainIndex)) %>% dplyr::rename (ID=`as.vector(trainIndex)`) 
    Test_ID  <<- cbind(Test,  as.vector(testIndex))  %>% dplyr::rename (ID=`as.vector(testIndex)`)  
    dplyr::all_equal(Train_ID, Test_ID)
    
    
    
    vars_forest <<- colnames( Train %>% dplyr::select(-c(1 , 2))  ) #extract varaible names by index, #excluded IV
    
    X <<- Train[,vars_forest]
    Y <<- as.vector( as.numeric( Train[,1] ) )
    W <<- as.vector( as.numeric( Train[,2] ) ) 
    #============================== raw full CF for training data ==============================
    cf_raw_key.tr <- CF_RAW_key(Train)   
    
    Y.hat  <<- cf_raw_key.tr$Y.hat
    W.hat  <<- cf_raw_key.tr$W.hat
    HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw
    varimp_cf  <- cf_raw_key.tr$varimp_cf
    selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx
    leafsize <<- LEAFSIZE_tune(round(Ntrain,0), treeNo, iterationNo)
    time_rawCF <- cf_raw_key.tr$time_rawCF
    
    
    ################Test perofrmance of iCF, oneCF on testing data #################
    #ATE, i.e.,no subgroup, just one overall group for all population
    ATE_all_tr <<- (DltY_DATA("NA",     Train_ID, "predict" ) %>%
                      dplyr::select(SubgroupID, Definition, dltY_ate) %>%
                      group_by(SubgroupID) %>%
                      summarize(dltY_ate = mean(dltY_ate)) )[1,1]
    ATE_all_te <<- (DltY_DATA("NA",     Test_ID, "predict" ) %>%
                      dplyr::select(SubgroupID, Definition, dltY_ate) %>%
                      group_by(SubgroupID) %>%
                      summarize(dltY_ate = mean(dltY_ate)) )[1,1]
    #---------------------------------------------------------------------------------------------------------------------
    # combine Delta Y results in a data frame for tree-based methods (explicitly show subgrup deicison, thus calcualte subgroup-specific treatment effect)
    #-------------------------------------------------------------------------------------------------------------------
    DltY_true_te  <<- DltY_DATA(tree_true_subgroup, Test_ID, "truth" )
    DltY_truth_te_max  <<- DltY_true_te %>%
      dplyr::select(SubgroupID, Definition, dltY_ate_t) %>% 
      dplyr::rename(dltY_ate = dltY_ate_t)
    CATE_max_true_te   <<- CATE_MAX_SG (DltY_truth_te_max )
    
    #vars_IV=c("X7")
    ###########################################
    # tune leaf size to grow D4, D3, D2 forests
    ###########################################
    #######################################################################################
    #try D2 first then times 2, 4, respectively to get the # for D3, D4 tree because:
    #1) it's simple: D2 tree has 1 split and 3 nodes; 2)if min.leaf.size too large, then D2 tree won't split, which leads to (aggregrate) error in the "vote_D2_tree", but D3, D4 tree won't lead to error
    #2) D4 sometimes has the same tree structure as D3, i.e. N of nodes = 7, thus should not start from D4 or count D4 N of nodes, only check D2 and D3.

# PS matched cohort for non-CF methods
#------------------------------------------
#traditional logistic regression for PS
#------------------------------------------
psform1 <<- as.formula(paste0("W ~ ", paste0(colnames( Train[,vars_forest] ), collapse = " + ")))
#read coefficient!!!, if a column has only one level (e.g. all 0 or all 1), the coefficient is "NA" and will follow the warning message when predict!
#Warning message:
#In predict.lm(object, newdata, se.fit, scale = 1, type = if (type ==  :
#                                                             prediction from a rank-deficient fit may be misleading
(aGlm_tr <- glm(psform1, family=binomial(), data=Train))
ps_tr <- predict(object = aGlm_tr, type = "response")
#------------------------------------------
#Logistic regression with lasso penalty doesn't work for next step, optmatch
#No need to 
#------------------------------------------
# fit propensity score model
#ps <- prop.func(as.matrix(X), W), 

ps_m_tr <- optmatch::match_on(aGlm_tr, caliper = 0.25*sd(ps_tr))
fm_tr   <- optmatch::fullmatch(ps_m_tr, data = Train)
Train_m <- na.omit(cbind(Train, fm_tr)) %>% select(-c(fm_tr))

table(Train_m$W)

(aGlm_te <- glm(psform1, family=binomial(), data=Test_ID))
ps_te <- predict(object = aGlm_te, type = "response")

ps_m_te <- optmatch::match_on(aGlm_te, caliper = 0.25*sd(ps_te))
fm_te   <- optmatch::fullmatch(ps_m_te, data = Test_ID)
#Test_ID_m <- na.omit(cbind(Test_ID, fm_te)) %>% select(-fm_te)




# ----------------------------------------------------
# SINGLE IT TREE ANALYSIS: only works for continuous outcome
# ----------------------------------------------------
s_IT =  Sys.time()
#__________________________
#if continuous outcome
if (length(unique(Y)) >=8){
  # OBTAIN A LARGE INITIAL TREE
  split.var <- 3:12  # COLUMNS OF COVARIATES
  col.y <- 1; col.trt <- 2; cols.split.var <- 3:12; 
  cols.nominal <- NULL
  # GROWING IT AND BOOTSTRAP PRUNING
  prn.boot <- bootstrap.grow.prune(B=50, 
                                   data=  Train_m %>%  dplyr::rename(trt=W, y=Y), 
                                   min.node.size=10, 
                                   n0=3, 
                                   max.depth=8,
                                   SSS=TRUE, 
                                   a=10, 
                                   q=0.01, 
                                   n.intervals=1,
                                   col.y=col.y, 
                                   col.trt=col.trt, 
                                   cols.split.var=cols.split.var, 
                                   cols.nominal=cols.nominal,
                                   mtry=10, 
                                   LeBlanc=TRUE, 
                                   min.boot.tree.size=1)  
  prn.result <- bootstrap.size(prn.boot, n=NROW(dat), plot.it= T)
  prn.result$G.a 
  # THE BEST IT TREE WITH a=2
  btree_it <- prn.result$btree[[1]] #lambda=2, which is AIC
  #plot.tree(btree_it, cols.nominal=7, textDepth=5, lines="rectangle")
  IT_btree_df <- CONVERT_tree_IT2CF(btree_it)
  #==================================================================================================================================================================
  if(nrow(IT_btree_df)==1){
    IT_subgroup="NA" 
    IT_subgroup_key= "NA" 
    vote_IT_subgroup= "NA" 
  } else {
    # IT_subgroup      <- PRE_MAJORITY_SUBGROUP ( list(IT_btree_df) );  
    #  IT_subgroup_key  <- lapply(IT_subgroup, function(df) subset(df, select=c("subgroupID", "subgroup")))
    #  vote_IT_subgroup <- MAJORITY_VOTE(list(IT_btree_df),  IT_subgroup_key,   IT_subgroup,   IT_subgroup, #as IT tree has not split frequency, this IT_subgroup is just to fool the Fx
    #                                    tree_true_subgroup, tree_true_N1,   tree_true_N123, split_val_round_posi)
    vote_IT_subgroup <- TREE2SUBGROUP(IT_btree_df)
  }
  #__________________________
  #if binary outcome
} else if (length(unique(Y)) == 2){
  IT_subgroup="NA" 
  IT_subgroup_key= "NA" 
  vote_IT_subgroup= "NA" 
} 

Deci_IT    <- OTHERTREE_DECI(vote_IT_subgroup)
SG2INT_IT <- SG2INT(Deci_IT)

e_IT =  Sys.time()
time_IT = difftime(e_IT, s_IT, units="secs")

disco_SG_IT     <- DISCOVERYRATE(Deci_IT,         truth_description, tree_true_subgroup, truth_INT, "IT", "subgroup")$value
disco_SG_IT.aa  <- DISCOVERYRATE(Deci_IT,         truth_description, tree_true_subgroup, truth_INT, "IT", "subgroup")$sg_Nacc_Nall
disco_SG_IT.aat <- DISCOVERYRATE(Deci_IT,         truth_description, tree_true_subgroup, truth_INT, "IT", "subgroup")$sg_Nacc_NallT
disco_INT_IT    <- DISCOVERYRATE(SG2INT(Deci_IT), truth_description, tree_true_subgroup, truth_INT, "IT", "interaction")$value

DltY_IT_te    <- DltY_DATA(Deci_IT,            Test_ID, "predict" )
disco_CATEmax_IT      <- DISCOVERYRATE(CATE_MAX_SG (DltY_IT_te),     truth_description, tree_true_subgroup, CATE_max_true_te, "IT", "CATEmax")$value
MSE_IT_te     <- MSE_DltY(DltY_true_te, DltY_IT_te)
MSE_ate_IT_te     <- as.numeric(MSE_IT_te$MSE_ate)
MSE_att_IT_te     <- as.numeric(MSE_IT_te$MSE_att)

############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm( prn.boot,   prn.result,  btree_it, IT_btree_df)
############################################################################
# ----------------------------------------------------
# aVirtualTwin ANALYSIS: only works for binary outcome
# ----------------------------------------------------
s_VT =  Sys.time()

if (length(unique(Y)) ==2){
  vt.o <- vt.data(Train_m, # %>% dplyr::mutate(Y= ifelse(Y==1,0,1)), 
                  "Y", 
                  "W", 
                  interactions = TRUE)
  set.seed(123)
  model.rf <- randomForest::randomForest(x = vt.o$getX(interactions = T),
                                         y = vt.o$getY(),
                                         ntree = 500)
  vt.f.rf <- vt.forest("one", vt.data = vt.o, 
                       model = model.rf, 
                       interactions = T)
  # grow RF for T = 1
  model.rf.trt1 <- randomForest(x = vt.o$getX(trt = 1), 
                                y = vt.o$getY(trt = 1))
  # grow RF for T = 0
  model.rf.trt0 <- randomForest(x = vt.o$getX(trt = 0), y = vt.o$getY(trt = 0))
  # initialize VT.forest.double()
  vt.doublef.rf <- vt.forest("double",
                             vt.data = vt.o, 
                             model_trt1 = model.rf.trt1, 
                             model_trt0 = model.rf.trt0)
  #model.fold <- vt.forest("fold", vt.data = vt.o, fold = 5, ratio = 1, interactions = T, ntree = 200)
  
  
  # initialize classification tree
  tr.class <- vt.tree("class",
                      #vt.difft = vt.f.rf,
                      vt.difft = vt.doublef.rf,
                      sens = ">",
                      threshold = quantile(vt.f.rf$difft, seq(.5, .8, .1)),
                      maxdepth = 3,
                      cp = 0,
                      maxcompete = 2) 
  # tr.class is a list if threshold is a vectoor
  class(tr.class)
  class(tr.class$tree1)
  
  tr.reg <- vt.tree("reg",
                    # vt.difft = vt.f.rf,
                    vt.difft = vt.doublef.rf,
                    sens = ">",
                    threshold = quantile(vt.f.rf$difft, seq(.5, .8, .1)))
  # tr.class is a list if threshold is a vectoor
  class(tr.reg)
  class(tr.reg$tree1)
  
  vt.sbgrps <- vt.subgroups(tr.class)
  
  # print tables with knitr package
  knitr::kable(vt.sbgrps)
  
  Deci_VT0 <- vt.sbgrps %>%
    dplyr::mutate(Subgroup_size = as.numeric(vt.sbgrps$`Subgroup size`)) %>%
    dplyr::filter(Subgroup_size > 50) %>% #require sample size > 50
    `rownames<-`( NULL ) %>% #remove original long rowname from VirtualTwin
    dplyr::select(Subgroup) %>%
    `colnames<-`( NULL ) 
  #modity the way to present binary (0,1) split value
  #remove "&" to prepare for sorting
  #dplyr::mutate(subgroup2 = stringr::str_replace_all(Subgroup, " & ", " "))  %>%
  #sort conditions in each subgroup definition
  #convert df to list 
  Deci_VT_L <-  split(Deci_VT0, seq(nrow(Deci_VT0)))
  #apply function that convert VirtualTwin format to Causal forest format
  Deci_VT_cf.format <- lapply(Deci_VT_L, VT2CF_format )
  
  Deci_VT_con_split2 <- data.frame(matrix(unlist(Deci_VT_cf.format), nrow=length(Deci_VT_cf.format), byrow=TRUE)) %>%
    `colnames<-`("Subgroup") %>% 
    #sort conditions in each subgroup definition
    rowwise() %>% 
    mutate(subgroup = paste(sort(unlist(strsplit(as.character(Subgroup), " & ", fixed = TRUE)), decreasing = TRUE), collapse = " & ")) 
  #sort subgroup definitions
  Deci_VT <-  Deci_VT_con_split2[order(Deci_VT_con_split2$subgroup),] %>%
    tibble::rowid_to_column() %>%
    dplyr::select(subgroup, rowid)%>%
    rename(subgroupID =rowid) %>%
    select(subgroupID, everything())
  
  
} else if (length(unique(Y)) >8){
  Deci_VT = "NA" 
} 
SG2INT_VT <- SG2INT(Deci_VT)

e_VT =  Sys.time()
time_VT = difftime(e_VT, s_VT, units="secs")

disco_SG_VT      <- DISCOVERYRATE(Deci_VT,         truth_description, tree_true_subgroup, truth_INT, "VT", "subgroup")$value
disco_SG_VT.aa   <- DISCOVERYRATE(Deci_VT,         truth_description, tree_true_subgroup, truth_INT, "VT", "subgroup")$sg_Nacc_Nall
disco_SG_VT.aat  <- DISCOVERYRATE(Deci_VT,         truth_description, tree_true_subgroup, truth_INT, "VT", "subgroup")$sg_Nacc_NallT
disco_INT_VT    <- DISCOVERYRATE(SG2INT(Deci_VT), truth_description, tree_true_subgroup, truth_INT, "VT", "interaction")$value

DltY_VT_te    <- DltY_DATA(Deci_VT,   Test_ID,   "predict" )
disco_CATEmax_VT      <- DISCOVERYRATE(CATE_MAX_SG (DltY_VT_te),     truth_description, tree_true_subgroup, CATE_max_true_te, "VT", "CATEmax")$value
MSE_VT_te     <- MSE_DltY(DltY_true_te, DltY_VT_te)
MSE_ate_VT_te     <- as.numeric(MSE_VT_te$MSE_ate)
MSE_att_VT_te     <- as.numeric(MSE_VT_te$MSE_att)


############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm( vt.o, model.rf, vt.f.rf, model.rf.trt1,  model.rf.trt0,  vt.doublef.rf,  tr.class, tr.reg, vt.sbgrps)
############################################################################


# --------------------------
# LASSO
# --------------------------
s_lasso =  Sys.time()
#not sure how to automatically use column name like (e.g., paste0( colnames(X), collapse = " + ")  ) to make this mmatrix, will figure out later:
mmatrix <- model.matrix(Y~(W + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)^4,Train_m)

colnames(mmatrix)
library(glmnet)
#if continuous outcome
if (length(unique(Y)) >=8){
  cvfit <- glmnet::cv.glmnet(x=mmatrix, y=Train_m$Y,alpha=1, family="gaussian")
  #if binary outcome
} else if (length(unique(Y)) ==2) {
  cvfit <- glmnet::cv.glmnet(x=mmatrix, y=Train_m$Y,alpha=1, family="binomial")
}
############################## Start Kevin's Edits #########################################

# Plot error as "function" of lambda
#plot(cvfit)
# Provides two "best" lambdas as dotted lines
# 1. lambda.min = lambda which minimizes CV mean squared error
# 2. lambda.1se = lambda with CV error within 1 SE of min AND has least nonzero parameters
# Next two lines plot the same thing, without transforming to log scale, marking only
# lambda.min
#plot(x=cvfit$lambda, y=cvfit$cvm)
#abline(v=cvfit$lambda[which(cvfit$cvm==min(cvfit$cvm))])

# Extract lambda with min cross validated mean error, look at coefficients under this lambda
lambda_min <- cvfit$lambda[which(cvfit$cvm==min(cvfit$cvm))]
which(cvfit$lambda == lambda_min)
coef(cvfit, s = lambda_min)
sum(coef(cvfit, s = lambda_min)!=0) # looks at number of non-zero parameters
# NOTE: If analyzing actual dataset, we don't know the truth, so this is the method I would 
# use to select lambda in a non-simulated dataset

lambda_min_results_pre <- data.frame("lambda"=lambda_min,
                                     "log_lambda"=log(lambda_min),
                                     "variable"=names(coef(cvfit, s = lambda_min)[,1]),
                                     "coef_ests"=coef(cvfit, s = lambda_min)[,1]) %>% 
  dplyr::filter(str_count(variable, ":") >=1 & str_count(variable, "W") ==1) %>%
  dplyr::mutate( coef_ests_abs=abs(coef_ests))   %>% 
  dplyr::arrange( desc(coef_ests_abs))  %>% 
  dplyr::mutate( se_coe = sd(coef_ests_abs)/sqrt(length(coef_ests_abs)),
                 mean_coe=mean(coef_ests_abs),
                 median_coe=median(coef_ests_abs),
                 max01_coe=max(coef_ests_abs)*0.5,
                 q1_coe=as.numeric(quantile(coef_ests_abs, 0.25)),
                 q3_coe=as.numeric(quantile(coef_ests_abs, 0.75)),
                 variable_int=as.character(variable)) 

lambda_min_results_0    <- lambda_min_results_pre %>% dplyr::filter (coef_ests_abs >  0)  
lambda_min_results_q1   <- lambda_min_results_pre %>% dplyr::filter (coef_ests_abs >  q1_coe)  
lambda_min_results_q2   <- lambda_min_results_pre %>% dplyr::filter (coef_ests_abs >  median_coe)
lambda_min_results_mean <- lambda_min_results_pre %>% dplyr::filter (coef_ests_abs >  mean_coe)  
lambda_min_results_q3   <- lambda_min_results_pre %>% dplyr::filter (coef_ests_abs >  q3_coe)  
lambda_min_results      <- lambda_min_results_pre %>% dplyr::filter (coef_ests_abs >  max01_coe)  


INT_lasso_0    <- EXTRACT_INT(lambda_min_results_0)
INT_lasso_q1   <- EXTRACT_INT(lambda_min_results_q1)
INT_lasso_q2   <- EXTRACT_INT(lambda_min_results_q2)
INT_lasso_mean <- EXTRACT_INT(lambda_min_results_mean)
INT_lasso_q3   <- EXTRACT_INT(lambda_min_results_q3)
INT_lasso      <- EXTRACT_INT(lambda_min_results)


e_lasso =  Sys.time()
time_lasso = difftime(e_lasso, s_lasso, units="secs")

disco_INT_lasso_0     <- DISCOVERYRATE(INT_lasso_0,    truth_description, tree_true_subgroup, truth_INT, "lasso",  "interaction")$value
disco_INT_lasso_q1    <- DISCOVERYRATE(INT_lasso_q1,   truth_description, tree_true_subgroup, truth_INT, "lasso",  "interaction")$value
disco_INT_lasso_q2    <- DISCOVERYRATE(INT_lasso_q2,   truth_description, tree_true_subgroup, truth_INT, "lasso",  "interaction")$value
disco_INT_lasso_mean  <- DISCOVERYRATE(INT_lasso_mean, truth_description, tree_true_subgroup, truth_INT, "lasso",  "interaction")$value
disco_INT_lasso_q3    <- DISCOVERYRATE(INT_lasso_q3,   truth_description, tree_true_subgroup, truth_INT, "lasso",  "interaction")$value
disco_INT_lasso       <- DISCOVERYRATE(INT_lasso,      truth_description, tree_true_subgroup, truth_INT, "lasso",  "interaction")$value

############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm(mmatrix, cvfit, lambda_min_results_pre)
############################################################################

# --------------------------
# FindIt (L2-SVM)
# --------------------------
s_FindIt =  Sys.time()

FindIt  <-FindIt( model.treat= Y ~ W,
                  model.main= ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                  model.int= ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                  data = Train_m,
                  type = ifelse(length(unique(Y)) >=8,  "continuous", "binary"), #cann't handle ordinal outcome with 3 levels, will worry about it later
                  treat.type="single")

FindIt_df <- as.data.frame( coef(summary(FindIt)) )

FindIt_results_pre <- data.frame(    variable  = rownames(FindIt_df),
                                     coef_ests = FindIt_df    ) %>% 
  dplyr::rename(coef_ests = 2) %>%
  dplyr::filter(str_count(variable, ":") >=1 & str_count(variable, "treat") ==1) %>%
  dplyr::mutate( coef_ests=as.numeric(coef_ests),
                 coef_ests_abs=abs(coef_ests))   %>% 
  dplyr::arrange ( desc(coef_ests_abs))  %>% 
  dplyr:: mutate( se_coe = sd(coef_ests_abs)/sqrt(length(coef_ests_abs)),
                  mean_coe=mean(coef_ests_abs),#10/11/2021, changed from "mean" to "median"
                  median_coe=median(coef_ests_abs),
                  max01_coe=max(coef_ests_abs)*0.5,
                  q1_coe=as.numeric(quantile(coef_ests_abs, 0.25)),
                  q3_coe=as.numeric(quantile(coef_ests_abs, 0.75)))  %>%
  tibble::remove_rownames() %>%
  dplyr::mutate (variable_c=as.character(variable)) %>%
  dplyr::mutate (variable_r = substring(variable_c, 6)) %>%
  dplyr::mutate (variable_int = paste0("W", variable_r))

FindIt_results_0    <- FindIt_results_pre %>% dplyr::filter (coef_ests_abs >   0)  
FindIt_results_q1   <- FindIt_results_pre %>% dplyr::filter (coef_ests_abs >   q1_coe)  
FindIt_results_q2   <- FindIt_results_pre %>% dplyr::filter (coef_ests_abs >   median_coe)  
FindIt_results_mean <- FindIt_results_pre %>% dplyr::filter (coef_ests_abs >   mean_coe)  
FindIt_results_q3   <- FindIt_results_pre %>% dplyr::filter (coef_ests_abs >   q3_coe)  
FindIt_results      <- FindIt_results_pre %>% dplyr::filter (coef_ests_abs >   max01_coe)  

#dplyr::filter (coef_ests_abs >  median_coe) %>%

#don't use ifelse, which does not return vectors!!!
INT_FindIt_0    <- EXTRACT_INT(FindIt_results_0)
INT_FindIt_q1   <- EXTRACT_INT(FindIt_results_q1)
INT_FindIt_q2   <- EXTRACT_INT(FindIt_results_q2)
INT_FindIt_mean <- EXTRACT_INT(FindIt_results_mean)
INT_FindIt_q3   <- EXTRACT_INT(FindIt_results_q3)
INT_FindIt      <- EXTRACT_INT(FindIt_results)

e_FindIt =  Sys.time()
time_FindIt = difftime(e_FindIt, s_FindIt, units="secs")

disco_INT_FindIt_0    <- DISCOVERYRATE(INT_FindIt_0,    truth_description, tree_true_subgroup, truth_INT, "FindIt", "interaction")$value
disco_INT_FindIt_q1   <- DISCOVERYRATE(INT_FindIt_q1,   truth_description, tree_true_subgroup, truth_INT, "FindIt", "interaction")$value
disco_INT_FindIt_q2   <- DISCOVERYRATE(INT_FindIt_q2,   truth_description, tree_true_subgroup, truth_INT, "FindIt", "interaction")$value
disco_INT_FindIt_mean <- DISCOVERYRATE(INT_FindIt_mean, truth_description, tree_true_subgroup, truth_INT, "FindIt", "interaction")$value
disco_INT_FindIt_q3   <- DISCOVERYRATE(INT_FindIt_q3,   truth_description, tree_true_subgroup, truth_INT, "FindIt", "interaction")$value
disco_INT_FindIt      <- DISCOVERYRATE(INT_FindIt,      truth_description, tree_true_subgroup, truth_INT, "FindIt", "interaction")$value



############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm(Train_m, FindIt, FindIt_df, FindIt_results_pre )
############################################################################
# ----------------------------------------------------
# Regularized Outcome weighted estimaton
# ----------------------------------------------------
#2) negative outcomes may cause ill-behaved outcome weighted estimation (Chen 2017. A general statistical framework for subgroup identification and comparative treatment scoring. Biometrics. doi:10.1111/biom.12676)
s_ROWSi=  Sys.time()
#make a new W_a vector for A learning, treatment=1, comparator=-1:
W_a = W
replace(W_a,  W_a==0, -1)
Hommel_P_threshold <<- 1e-9
#=============================
sg_mod_ROWSi<- fit.subgroup(x = as.matrix(X), y = Y , trt = W_a,
                            propensity.func = prop.func,
                            method= "weighting", 
                            #"a_learning",
                            loss = #"owl_hinge_loss", super time-consuming for large dataset, thus give up!!!
                              #"sq_loss_lasso",
                              #"owl_logistic_loss_lasso",#ROWSi
                              "owl_logistic_flip_loss_lasso",#ROWSi, "flip" allow for non-positive outcomes, All losses for continuous outcomes can be used for binary outcome!!!
                            nfolds = 10)
#summary(sg_mod_owe)
#when 2W INT, error:
#Error in summarize.subgroups.default(x = xx, subgroup = subgroup) : 
#  There is only one unique subgroup. No subgroups to compare with.
if ( is.null( tryCatch(summarize.subgroups(sg_mod_ROWSi), error=function(e){}) ) == TRUE ) {
  cov_int_ROWSi_df="NA"
  INT_ROWSi_0="NA"
  INT_ROWSi_q1="NA"
  INT_ROWSi_q2="NA"
  INT_ROWSi_mean="NA"
  INT_ROWSi="NA"
  INT_ROWSi_q3="NA"
} else {
  comp_ROWSi<- summarize.subgroups(sg_mod_ROWSi)
  #given threshold of a multiple comparisons-adjsuted P-value
  #print only the covariates which have "significant' differences between subgroups with a multiple comparisons-adjusted p-value
  #-----------------------------------------------------------------
  #10/04/2021: previous way to identify interaction is not accurate, which actually identify potential subgroups defined by covariates by recommended therapy
  cov_int_ROWSi_df <- as.data.frame(print(comp_ROWSi, p.value = Hommel_P_threshold)) 
  
  cov_int_ROWSi_df2 <- cov_int_ROWSi_df %>%
    dplyr::mutate (variable_c = rownames(cov_int_ROWSi_df)) %>%
    dplyr::mutate (variable_int = paste0("W:", variable_c))
  #-----------------------------------------------------------------
  
  
  summ <- summary( sg_mod_ROWSi$coefficients)
  cov_int_ROWSi_df3 <- data.frame(variable     = rownames(sg_mod_ROWSi$coefficients)[summ$i],
                                  Destination = colnames(sg_mod_ROWSi$coefficients)[summ$j],
                                  coef_ests      = summ$x) %>%
    dplyr::filter(variable != "Trt1") %>% 
    dplyr::mutate( coef_ests_abs=abs(coef_ests))   %>% 
    dplyr::arrange ( desc(coef_ests_abs)) %>%
    dplyr::mutate (variable_int = paste0("W:", variable), 
                   coef_ests_r = round(coef_ests,3),
                   se_coe = sd(coef_ests_abs)/sqrt(length(coef_ests_abs)),
                   mean_coe=mean(coef_ests_abs),#10/11/2021, changed from "mean" to "median"
                   median_coe=median(coef_ests_abs),
                   max01_coe=max(coef_ests_abs)*0.5,
                   q1_coe=as.numeric(quantile(coef_ests_abs, 0.25)),
                   q3_coe=as.numeric(quantile(coef_ests_abs, 0.75))) 
  
  ROWSi_results_0    <- cov_int_ROWSi_df3 %>% dplyr::filter (coef_ests_abs >   0)  
  ROWSi_results_q1   <- cov_int_ROWSi_df3 %>% dplyr::filter (coef_ests_abs >   q1_coe)  
  ROWSi_results_q2   <- cov_int_ROWSi_df3 %>% dplyr::filter (coef_ests_abs >   median_coe)  
  ROWSi_results_mean <- cov_int_ROWSi_df3 %>% dplyr::filter (coef_ests_abs >   mean_coe)  
  ROWSi_results_q3   <- cov_int_ROWSi_df3 %>% dplyr::filter (coef_ests_abs >   q3_coe)  
  ROWSi_results      <- cov_int_ROWSi_df3 %>% dplyr::filter (coef_ests_abs >   max01_coe)  
  
  
  #don't use ifelse, which does not return vectors!!!
  INT_ROWSi_0    <- EXTRACT_INT(ROWSi_results_0 ) 
  INT_ROWSi_q1   <- EXTRACT_INT(ROWSi_results_q1)  
  INT_ROWSi_q2   <- EXTRACT_INT(ROWSi_results_q2 )
  INT_ROWSi_mean <- EXTRACT_INT(ROWSi_results_mean )  
  INT_ROWSi_q3   <- EXTRACT_INT(ROWSi_results_q3 )  
  INT_ROWSi      <- EXTRACT_INT(ROWSi_results )  
  
  
  
}
e_ROWSi=  Sys.time()
time_ROWSi=  difftime(e_ROWSi, s_ROWSi, units="secs")


disco_INT_ROWSi_0     <- DISCOVERYRATE(INT_ROWSi_0,     truth_description, tree_true_subgroup, truth_INT, "ROWSi",  "interaction")$value
disco_INT_ROWSi_q1    <- DISCOVERYRATE(INT_ROWSi_q1,    truth_description, tree_true_subgroup, truth_INT, "ROWSi",  "interaction")$value
disco_INT_ROWSi_q2    <- DISCOVERYRATE(INT_ROWSi_q2,    truth_description, tree_true_subgroup, truth_INT, "ROWSi",  "interaction")$value
disco_INT_ROWSi_mean  <- DISCOVERYRATE(INT_ROWSi_mean,  truth_description, tree_true_subgroup, truth_INT, "ROWSi",  "interaction")$value
disco_INT_ROWSi_q3    <- DISCOVERYRATE(INT_ROWSi_q3,    truth_description, tree_true_subgroup, truth_INT, "ROWSi",  "interaction")$value
disco_INT_ROWSi       <- DISCOVERYRATE(INT_ROWSi,       truth_description, tree_true_subgroup, truth_INT, "ROWSi",  "interaction")$value
######################################

############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm(W_a, comp_ROWSi, cov_int_ROWSi_df, cov_int_ROWSi_df2, cov_int_ROWSi_df3)
############################################################################





#%>%
# rename_at(vars(-SubID_iCF), ~ paste0(., '.ori')) %>% 
#rename(SubID_iCF.ori = SubID_iCF
############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm(Dlt_Y_iCF_AIC, SG_D_iCF_AIC)
############################################################################

#---------------------------------------------------------------------------------------------
#       oneCF CV
#---------------------------------------------------------------------------------------------
oneCFCV_SG <- iCFCV(10, treeNo, 1)
Deci_CV_oneCF.act.bias =oneCFCV_SG$seletedSG_CV.act.bias
Deci_CV_oneCF.tran.bias =oneCFCV_SG$seletedSG_CV.tran.bias
Deci_CV_oneCF.act.AIC =oneCFCV_SG$seletedSG_CV.act.AIC
Deci_CV_oneCF.tran.AIC =oneCFCV_SG$seletedSG_CV.tran.AIC
CVL_Deci_oneCF.act.AIC =  oneCFCV_SG$Deci_AIC_test.act
CVL_Deci_oneCF.tran.AIC = oneCFCV_SG$Deci_AIC_test.tran
time_oneCF_CV = oneCFCV_SG$time_iCFCV

###################
iterationNo
treeNo
Deci_CV_oneCF.act.bias
Deci_CV_oneCF.tran.bias
Deci_CV_oneCF.act.AIC
Deci_CV_oneCF.tran.AIC
time_oneCF_CV
CVL_Deci_oneCF.act.AIC
CVL_Deci_oneCF.tran.AIC
###################



if (Deci_CV_oneCF.act.bias=="NA"){
  oneCF_CV_per.act.bias  <- TREE_PERFORMANCE(Deci_CV_oneCF.act.bias,  1)  
} else {
  oneCF_CV_per.act.bias  <- TREE_PERFORMANCE(Deci_CV_oneCF.act.bias$majority,  1)
}

if (Deci_CV_oneCF.tran.bias=="NA"){
  oneCF_CV_per.tran.bias <- TREE_PERFORMANCE(Deci_CV_oneCF.tran.bias, 1)
} else {
  oneCF_CV_per.tran.bias <- TREE_PERFORMANCE(Deci_CV_oneCF.tran.bias$majority, 1)
}

oneCF_CV_per.act.AIC  <- TREE_PERFORMANCE(Deci_CV_oneCF.act.AIC,  1)
oneCF_CV_per.tran.AIC <- TREE_PERFORMANCE(Deci_CV_oneCF.tran.AIC, 1)

disco_SG_CV_oneCF.act.bias       =  oneCF_CV_per.act.bias$disco_SG_iCF 
disco_INT_CV_oneCF.act.bias      =  oneCF_CV_per.act.bias$disco_INT_iCF 
disco_CATEmax_CV_oneCF.act.bias  =  oneCF_CV_per.act.bias$disco_CATEmax_iCF
MSE_ate_CV_oneCF_te.act.bias     =  oneCF_CV_per.act.bias$MSE_ate_iCF_te
MSE_att_CV_oneCF_te.act.bias     =  oneCF_CV_per.act.bias$MSE_att_iCF_te

disco_SG_CV_oneCF.tran.bias       =  oneCF_CV_per.tran.bias$disco_SG_iCF 
disco_INT_CV_oneCF.tran.bias      =  oneCF_CV_per.tran.bias$disco_INT_iCF 
disco_CATEmax_CV_oneCF.tran.bias  =  oneCF_CV_per.tran.bias$disco_CATEmax_iCF
MSE_ate_CV_oneCF_te.tran.bias     =  oneCF_CV_per.tran.bias$MSE_ate_iCF_te
MSE_att_CV_oneCF_te.tran.bias     =  oneCF_CV_per.tran.bias$MSE_att_iCF_te

disco_SG_CV_oneCF.act.AIC       =  oneCF_CV_per.act.AIC$disco_SG_iCF 
disco_INT_CV_oneCF.act.AIC      =  oneCF_CV_per.act.AIC$disco_INT_iCF 
disco_CATEmax_CV_oneCF.act.AIC  =  oneCF_CV_per.act.AIC$disco_CATEmax_iCF
MSE_ate_CV_oneCF_te.act.AIC     =  oneCF_CV_per.act.AIC$MSE_ate_iCF_te
MSE_att_CV_oneCF_te.act.AIC     =  oneCF_CV_per.act.AIC$MSE_att_iCF_te

disco_SG_CV_oneCF.tran.AIC       =  oneCF_CV_per.tran.AIC$disco_SG_iCF 
disco_INT_CV_oneCF.tran.AIC      =  oneCF_CV_per.tran.AIC$disco_INT_iCF 
disco_CATEmax_CV_oneCF.tran.AIC  =  oneCF_CV_per.tran.AIC$disco_CATEmax_iCF
MSE_ate_CV_oneCF_te.tran.AIC     =  oneCF_CV_per.tran.AIC$MSE_ate_iCF_te
MSE_att_CV_oneCF_te.tran.AIC     =  oneCF_CV_per.tran.AIC$MSE_att_iCF_te

#=========================================================================================================
#                                            oneCF by AIC
#=========================================================================================================
Deci_NCV_oneCF.act.AIC  = oneCFCV_SG$seletedSG_NCV.act.AIC
Deci_NCV_oneCF.tran.AIC = oneCFCV_SG$seletedSG_NCV.tran.AIC

time_oneCF_NCV = oneCFCV_SG$time_iCFCV/10

oneCF_NCV_per.act.AIC  <- TREE_PERFORMANCE(Deci_NCV_oneCF.act.AIC,  1)
oneCF_NCV_per.tran.AIC <- TREE_PERFORMANCE(Deci_NCV_oneCF.tran.AIC, 1)

disco_SG_NCV_oneCF.act.AIC       =  oneCF_NCV_per.act.AIC$disco_SG_iCF 
disco_INT_NCV_oneCF.act.AIC      =  oneCF_NCV_per.act.AIC$disco_INT_iCF 
disco_CATEmax_NCV_oneCF.act.AIC  =  oneCF_NCV_per.act.AIC$disco_CATEmax_iCF
MSE_ate_NCV_oneCF_te.act.AIC     =  oneCF_NCV_per.act.AIC$MSE_ate_iCF_te
MSE_att_NCV_oneCF_te.act.AIC     =  oneCF_NCV_per.act.AIC$MSE_att_iCF_te

disco_SG_NCV_oneCF.tran.AIC       =  oneCF_NCV_per.tran.AIC$disco_SG_iCF 
disco_INT_NCV_oneCF.tran.AIC      =  oneCF_NCV_per.tran.AIC$disco_INT_iCF 
disco_CATEmax_NCV_oneCF.tran.AIC  =  oneCF_NCV_per.tran.AIC$disco_CATEmax_iCF
MSE_ate_NCV_oneCF_te.tran.AIC     =  oneCF_NCV_per.tran.AIC$MSE_ate_iCF_te
MSE_att_NCV_oneCF_te.tran.AIC     =  oneCF_NCV_per.tran.AIC$MSE_att_iCF_te

############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm(oneCFCV_SG)
############################################################################

#---------------------------------------------------------------------------------------------
#       iCF CV
#---------------------------------------------------------------------------------------------
iCFCV_SG <- iCFCV(10, treeNo, iterationNo, 0.1)





Deci_CV_iCF.act.bias = iCFCV_SG$seletedSG_CV.act.bias
Deci_CV_iCF.tran.bias =iCFCV_SG$seletedSG_CV.tran.bias
Deci_CV_iCF.act.AIC = iCFCV_SG$seletedSG_CV.act.AIC
Deci_CV_iCF.tran.AIC =iCFCV_SG$seletedSG_CV.tran.AIC
CVL_Deci_iCF.act.AIC =  iCFCV_SG$CVL_Deci_AIC_test.act
CVL_Deci_iCF.tran.AIC = iCFCV_SG$CVL_Deci_AIC_test.tran
V_D4_subgroup_iCFCV = iCFCV_SG$vote_D4_subgroup.L
V_D3_subgroup_iCFCV = iCFCV_SG$vote_D3_subgroup.L
V_D2_subgroup_iCFCV = iCFCV_SG$vote_D2_subgroup.L
leafsize_iCFCV = iCFCV_SG$leafsize
time_iCF_CV = iCFCV_SG$time_iCFCV

###################
iterationNo
treeNo
Deci_CV_iCF.act.bias
Deci_CV_iCF.tran.bias
Deci_CV_iCF.act.AIC
Deci_CV_iCF.tran.AIC
time_iCF_CV
CVL_Deci_iCF.act.AIC
CVL_Deci_iCF.tran.AIC
###################


if (Deci_CV_iCF.act.bias=="NA"){
iCF_CV_per.act.bias    <- TREE_PERFORMANCE(Deci_CV_iCF.act.bias,  iterationNo )
} else {
iCF_CV_per.act.bias    <- TREE_PERFORMANCE(Deci_CV_iCF.act.bias$majority,  iterationNo )
}

if (Deci_CV_iCF.tran.bias=="NA"){
iCF_CV_per.tran.bias   <- TREE_PERFORMANCE(Deci_CV_iCF.tran.bias, iterationNo )
} else {
iCF_CV_per.tran.bias   <- TREE_PERFORMANCE(Deci_CV_iCF.tran.bias$majority, iterationNo )
}

iCF_CV_per.act.AIC    <- TREE_PERFORMANCE(Deci_CV_iCF.act.AIC,  iterationNo )
iCF_CV_per.tran.AIC   <- TREE_PERFORMANCE(Deci_CV_iCF.tran.AIC, iterationNo )

disco_SG_CV_iCF.act.bias       =  iCF_CV_per.act.bias$disco_SG_iCF 
disco_INT_CV_iCF.act.bias      =  iCF_CV_per.act.bias$disco_INT_iCF 
disco_CATEmax_CV_iCF.act.bias  =  iCF_CV_per.act.bias$disco_CATEmax_iCF
MSE_ate_CV_iCF_te.act.bias     =  iCF_CV_per.act.bias$MSE_ate_iCF_te
MSE_att_CV_iCF_te.act.bias     =  iCF_CV_per.act.bias$MSE_att_iCF_te

disco_SG_CV_iCF.tran.bias       =  iCF_CV_per.tran.bias$disco_SG_iCF 
disco_INT_CV_iCF.tran.bias      =  iCF_CV_per.tran.bias$disco_INT_iCF 
disco_CATEmax_CV_iCF.tran.bias  =  iCF_CV_per.tran.bias$disco_CATEmax_iCF
MSE_ate_CV_iCF_te.tran.bias     =  iCF_CV_per.tran.bias$MSE_ate_iCF_te
MSE_att_CV_iCF_te.tran.bias     =  iCF_CV_per.tran.bias$MSE_att_iCF_te

disco_SG_CV_iCF.act.AIC       =  iCF_CV_per.act.AIC$disco_SG_iCF 
disco_INT_CV_iCF.act.AIC      =  iCF_CV_per.act.AIC$disco_INT_iCF 
disco_CATEmax_CV_iCF.act.AIC  =  iCF_CV_per.act.AIC$disco_CATEmax_iCF
MSE_ate_CV_iCF_te.act.AIC     =  iCF_CV_per.act.AIC$MSE_ate_iCF_te
MSE_att_CV_iCF_te.act.AIC     =  iCF_CV_per.act.AIC$MSE_att_iCF_te

disco_SG_CV_iCF.tran.AIC       =  iCF_CV_per.tran.AIC$disco_SG_iCF 
disco_INT_CV_iCF.tran.AIC      =  iCF_CV_per.tran.AIC$disco_INT_iCF 
disco_CATEmax_CV_iCF.tran.AIC  =  iCF_CV_per.tran.AIC$disco_CATEmax_iCF
MSE_ate_CV_iCF_te.tran.AIC     =  iCF_CV_per.tran.AIC$MSE_ate_iCF_te
MSE_att_CV_iCF_te.tran.AIC     =  iCF_CV_per.tran.AIC$MSE_att_iCF_te

#=========================================================================================================
#                                              iCF by AIC
#=========================================================================================================

Deci_NCV_iCF.act.AIC = iCFCV_SG$seletedSG_NCV.act.AIC
Deci_NCV_iCF.tran.AIC =iCFCV_SG$seletedSG_NCV.tran.AIC

time_iCF_NCV = iCFCV_SG$time_iCFCV/10

iCF_NCV_per.act.AIC    <- TREE_PERFORMANCE(Deci_NCV_iCF.act.AIC,  iterationNo )
iCF_NCV_per.tran.AIC   <- TREE_PERFORMANCE(Deci_NCV_iCF.tran.AIC, iterationNo )

disco_SG_NCV_iCF.act.AIC       =  iCF_NCV_per.act.AIC$disco_SG_iCF 
disco_INT_NCV_iCF.act.AIC      =  iCF_NCV_per.act.AIC$disco_INT_iCF 
disco_CATEmax_NCV_iCF.act.AIC  =  iCF_NCV_per.act.AIC$disco_CATEmax_iCF
MSE_ate_NCV_iCF_te.act.AIC     =  iCF_NCV_per.act.AIC$MSE_ate_iCF_te
MSE_att_NCV_iCF_te.act.AIC     =  iCF_NCV_per.act.AIC$MSE_att_iCF_te

disco_SG_NCV_iCF.tran.AIC       =  iCF_NCV_per.tran.AIC$disco_SG_iCF 
disco_INT_NCV_iCF.tran.AIC      =  iCF_NCV_per.tran.AIC$disco_INT_iCF 
disco_CATEmax_NCV_iCF.tran.AIC  =  iCF_NCV_per.tran.AIC$disco_CATEmax_iCF
MSE_ate_NCV_iCF_te.tran.AIC     =  iCF_NCV_per.tran.AIC$MSE_ate_iCF_te
MSE_att_NCV_iCF_te.tran.AIC     =  iCF_NCV_per.tran.AIC$MSE_att_iCF_te
############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
rm(iCFCV_SG)
############################################################################




#2nd combine with all other vectors: 
finaldata <- #cbind(
              as.data.frame( cbind(intTRUE, UC_indi, scenario,   j,  outcome_type, nstudy,
                                        time_iCF_CV ,
                                        time_oneCF_CV,
                                        time_iCF_NCV ,
                                        time_oneCF_NCV,

                                        disco_SG_CV_iCF.act.bias, 
                                        disco_INT_CV_iCF.act.bias,
                                        disco_CATEmax_CV_iCF.act.bias,
                                        MSE_ate_CV_iCF_te.act.bias,
                                        MSE_att_CV_iCF_te.act.bias,
                                        
                                        disco_SG_CV_iCF.tran.bias,
                                        disco_INT_CV_iCF.tran.bias,
                                        disco_CATEmax_CV_iCF.tran.bias,
                                        MSE_ate_CV_iCF_te.tran.bias,
                                        MSE_att_CV_iCF_te.tran.bias,
                                        
                                        disco_SG_CV_iCF.act.AIC,
                                        disco_INT_CV_iCF.act.AIC,
                                        disco_CATEmax_CV_iCF.act.AIC,
                                        MSE_ate_CV_iCF_te.act.AIC,
                                        MSE_att_CV_iCF_te.act.AIC,
                                        
                                        disco_SG_CV_iCF.tran.AIC,
                                        disco_INT_CV_iCF.tran.AIC,
                                        disco_CATEmax_CV_iCF.tran.AIC,
                                        MSE_ate_CV_iCF_te.tran.AIC,
                                        MSE_att_CV_iCF_te.tran.AIC,
                                   
                                       disco_SG_NCV_iCF.act.AIC,
                                       disco_INT_NCV_iCF.act.AIC,
                                       disco_CATEmax_NCV_iCF.act.AIC,
                                       MSE_ate_NCV_iCF_te.act.AIC,
                                       MSE_att_NCV_iCF_te.act.AIC,
                                       
                                       disco_SG_NCV_iCF.tran.AIC,
                                       disco_INT_NCV_iCF.tran.AIC,
                                       disco_CATEmax_NCV_iCF.tran.AIC,
                                       MSE_ate_NCV_iCF_te.tran.AIC,
                                       MSE_att_NCV_iCF_te.tran.AIC,

                                        disco_SG_CV_oneCF.act.bias,
                                        disco_INT_CV_oneCF.act.bias,
                                        disco_CATEmax_CV_oneCF.act.bias,
                                        MSE_ate_CV_oneCF_te.act.bias,
                                        MSE_att_CV_oneCF_te.act.bias,
                                        
                                        disco_SG_CV_oneCF.tran.bias,
                                        disco_INT_CV_oneCF.tran.bias ,
                                        disco_CATEmax_CV_oneCF.tran.bias,
                                        MSE_ate_CV_oneCF_te.tran.bias ,
                                        MSE_att_CV_oneCF_te.tran.bias,
                                        
                                        disco_SG_CV_oneCF.act.AIC,
                                        disco_INT_CV_oneCF.act.AIC,
                                        disco_CATEmax_CV_oneCF.act.AIC,
                                        MSE_ate_CV_oneCF_te.act.AIC,
                                        MSE_att_CV_oneCF_te.act.AIC,
                                        
                                        disco_SG_CV_oneCF.tran.AIC,
                                        disco_INT_CV_oneCF.tran.AIC,
                                        disco_CATEmax_CV_oneCF.tran.AIC,
                                        MSE_ate_CV_oneCF_te.tran.AIC,
                                        MSE_att_CV_oneCF_te.tran.AIC,
                                   
                                       disco_SG_NCV_oneCF.act.AIC,
                                       disco_INT_NCV_oneCF.act.AIC,
                                       disco_CATEmax_NCV_oneCF.act.AIC,
                                       MSE_ate_NCV_oneCF_te.act.AIC,
                                       MSE_att_NCV_oneCF_te.act.AIC,
                                       
                                       disco_SG_NCV_oneCF.tran.AIC,
                                       disco_INT_NCV_oneCF.tran.AIC,
                                       disco_CATEmax_NCV_oneCF.tran.AIC,
                                       MSE_ate_NCV_oneCF_te.tran.AIC,
                                       MSE_att_NCV_oneCF_te.tran.AIC,
                                        #------------------------------------------------------------
                                        # SUBGROUP demonstration
                                        #-----------------------------------------------------------
                                        #1st: HTE presense test P-value
                                        HTE_P_cf.raw, #raw CF P-value
                                        #2st: split frequency: make sense to get data across simulations? as when data changes, splitting variable changes
                                        
                                        #3rd: Delta Y in 2nd step
                                        time_rawCF,      
                    
                                        time_IT,      time_VT,   time_lasso,      time_FindIt ,      time_ROWSi, 

                                        disco_INT_IT   ,
                                        disco_INT_VT   ,
                                        
                                        disco_SG_IT, 
                                        disco_SG_VT, 
                                        
                                        disco_CATEmax_IT   ,
                                        disco_CATEmax_VT   ,
                                        
                                        disco_INT_lasso_0,  disco_INT_lasso_q1,  disco_INT_lasso_q2,  disco_INT_lasso_mean,  disco_INT_lasso_q3,  disco_INT_lasso,
                                        disco_INT_FindIt_0, disco_INT_FindIt_q1, disco_INT_FindIt_q2, disco_INT_FindIt_mean, disco_INT_FindIt_q3, disco_INT_FindIt,
                                        disco_INT_ROWSi_0,  disco_INT_ROWSi_q1,  disco_INT_ROWSi_q2,  disco_INT_ROWSi_mean,  disco_INT_ROWSi_q3,  disco_INT_ROWSi,

                                        MSE_ate_IT_te,
                                        MSE_att_IT_te,
                                        
                                        MSE_ate_VT_te,
                                        MSE_att_VT_te
                                        
))#, 


save(finaldata, file=#paste0("finaldata_", rep, ".RData")
     
                    paste0( paste(UC_indi,  outcome_type, treeNo, iterationNo, nstudy, intTRUE,  sep="_") , rep,  ".RData" ) )




#------------------------------------------------------------------
# 2nd step: combine again with dataframe
#------------------------------------------------------------------
#Dlt_Y_iCF_AIC_df2 
#)







##################################################################################################
#############################################Combining Data AND Writing Data to File##############
##################################################################################################
# get command line arguments
#args = commandArgs(trailingOnly = TRUE)
#rep = as.numeric(args[1])

#set.seed(20160413 + rep)
#save(finaldata, file="asthma_glpvdpp_2y_r.RData")


#if ( any(duplicated(names(finaldata)))  == TRUE ) {
#  warning("Duplicate column name of finaldata, quit R!")
#  q()
#} else {
  
#  if(j==1){
#    write.table(finaldata,file=paste(UC_indi,  outcome_type, treeNo, iterationNo, nstudy, intTRUE, ".csv", sep="_"), row.names=FALSE,col.names=TRUE,    sep=",",append=TRUE)
#  } else{
#    write.table(finaldata,file=paste(UC_indi,  outcome_type, treeNo, iterationNo, nstudy, intTRUE, ".csv", sep="_"), row.names=FALSE,col.names=FALSE,   sep=",",append=TRUE)
#  }
#  
#}

############################################################################
#        RELEASE MEMORY I : remove big dataset/vector to speed up R        #
#rm(Dlt_Y_iCF_AIC_df2)
############################################################################



	} #this } for loop
  
  # sink()
  
}   #this } for sim.function

#sim.function(scenario,	intTRUE, UCindi, outcome_type,   a8,	 a9,	a10,  a11,  a12, a12_, a13, a13_,  a14, a14_ ,a15, a15_,  a16, a16_ , a17, a17_, a18, a18_, a21, a22, a23, a2s2, a4s2, a1i3,  a2i4,  a2i3,a3i4,  treeNo, iterationNo, nstart, nsim, nstudy, split_val_round_posi)
sim.function('A', '0W_000_00_00',"NoUCCV", "continuous", 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,      0,    0,   0,   0,    0,   0,     0,   0,    0,   0,   0,    0,   0,     0,   0,    0,       200,     50,        1,     1,    10000, 1) 
