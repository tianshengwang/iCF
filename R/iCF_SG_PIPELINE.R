
###############################################################################################
# causal forest wrapper function for easy-input of different parameter
###############################################################################################
#' CF
#' 
#' This function is the wrapper function to run raw causal forest with predefined min.node.size = round(nrow(Train)/depth, 0), 
#' @param depth parameter to tune min.node.size = round(nrow(Train)/depth, 0), 
#' @param treeNo number of trees in each forest
#' @param tunepara request wihout tuning
#' 
#' @return the key results from raw causal forest 
#'
#' @export
#' 
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
#' iCF
#' 
#' function to run iterative causal forest at a certain depth
#' @param depth parameter to tune min.node.size = round(nrow(Train)/depth, 0), 
#' @param treeNo number of trees in each forest
#' @param tunepara request wihout tuning
#' @param iterationNo iteration No (if = 1 then oneCF)
#' @param tree_depth "D2", "D3", "D4", or "D5"
#' @param split_val_round_posi split value rounding position
#'  
#' @return final subgroup decision G_iCF
#' 
#' @export
#' 
iCF <- function( depth, treeNo, tunepara = "none", iterationNo, tree_depth, split_val_round_posi){
  besttreelist = list() #making an empty list
  besttreelist_L = list () #for list format (rather than df format) of best trees
  splitfreqlist = list()
  treeBlist <- list()
  #treeBlist_L <- list()
  mvtreelist = list()
  # mvtreelist_L = list()
  for (k in 1:iterationNo) { 
    cf <- CF(depth, # =as.numeric(sample(leafsize, 1)), 
             treeNo, 
             tunepara)
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
#' oneCF
#' 
#' function to run iterative causal forest at a certain depth
#' @param depth parameter to tune min.node.size = round(nrow(Train)/depth, 0), 
#' @param treeNo number of trees in each forest
#' @param tunepara request wihout tuning
#' @param iterationNo iteration No (if = 1 then oneCF)
#' @param tree_depth "D2", "D3", "D4", or "D5"
#' @param split_val_round_posi split value rounding position
#'  
#' @return final subgroup decision G_iCF
#' 
#' @export

oneCF <- function(depth, treeNo, tunepara="none", tree_depth, split_val_round_posi){
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


#' find_level_over2
#' 
#' function to find variables with > 2 levls
#' @param X features
#'  
#' @return varialbes more than 2 levels in X
#' 
#' @export

find_level_over2 <- function(X_dat){
vars_levelover2 <-sapply(X_dat, function(x) length(unique(x)) )%>% as.data.frame() %>% 
  dplyr::rename(level = 1) %>% tibble::rownames_to_column() %>%
  dplyr::filter(level>2 & level <=8) 

if (nrow(vars_levelover2)==0) {
  var_catover2 <<- NA 
} else {
  var_catover2 <<-vars_levelover2$rowname
}

return(var_catover2)
}



#' SUBGROUP_PIPELINE
#' 
#' function to get subroup decision G_D, get ready for predicting Y (transformed outcome) in testing set,  and models build from G_D to predict Y*
#' @param X features
#' @param Y outcome
#' @param W exposure 
#' @param Y.hat predicted outcome
#' @param W.hat predicted propensity score 
#' @param treeNo tree No
#' @param iterationNo iteration No
#' @param Ntrain sample size of training set for developing G_D
#' @param Train_ID training set in Cross-Validation with the CV ID
#' @param Test_ID testing set in Cross-Validation with the CV ID
#' @param variable_type high-dimensional (HD) or non-HD
#' @param HTE_P_cf.raw P-value for HTE test in raw CF
#' @param P_threshold threshold for P value
#' @param split_val_round_posi split value rounding position
#'  
#' @return G_D
#' 
#' @export
 
SUBGROUP_PIPELINE<- function(X, 
                            Y, 
                            W, 
                            Y.hat,
                            W.hat,
                            selected_cf.idx,
                            leafsize, 
                            treeNo, 
                            iterationNo, 
                            Ntrain, 
                            Train_ID, 
                            Test_ID, 
                            variable_type,
                            HTE_P_cf.raw ,
                            P_threshold, 
                            split_val_round_posi){
  
  s_iCF=Sys.time()
  
  if (iterationNo>1){
    #Green Column in iCF Algorithm Figure
    #when leaf size too small, D2 has this error "Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) : no rows to aggregate"
    #thus to make it auto run, if iCF_D2 has error, then make it equal to iCF_D3 temparily.
    #after tuning leaf size, this erorr unlikely to happen

    iCF_D5<- iCF(leafsize$D5,  treeNo, tunepara = "none", iterationNo,  "D5",  split_val_round_posi)
    
    iCF_D4<- iCF(leafsize$D4,  treeNo, tunepara = "none", iterationNo,  "D4",  split_val_round_posi)
    
    iCF_D3<- iCF(leafsize$D3,  treeNo, tunepara = "none", iterationNo,  "D3",  split_val_round_posi) 
    
    if ( is.null( tryCatch( iCF(leafsize$D2,  treeNo, tunepara = "none", iterationNo,  "D2",  split_val_round_posi), error=function(e){}) ) == TRUE ) {
      warning("The minimum leaf size is too small for Depth 2 causal forest, need to increase it!")
      iCF_D2 = iCF_D3
    } else {
      iCF_D2 = iCF(leafsize$D2,  treeNo, tunepara = "none", iterationNo, "D2",  split_val_round_posi)
    }
   
    
    #best tree dataframe, Blue Column in iCF Algorithm Figure
    iCF_D5_BT <<- GET_TREE_L(iCF_D5, 1) 
    iCF_D4_BT <<- GET_TREE_L(iCF_D4, 1) 
    iCF_D3_BT <<- GET_TREE_L(iCF_D3, 1)
    iCF_D2_BT <<- GET_TREE_L(iCF_D2, 1) 
    
    #mean number of node from a causal tree dataframe
    mean_node_D5 <-round(mean( do.call("rbind" , lapply(iCF_D5_BT, function(df) nrow(df))) ), 0)
    mean_node_D4 <-round(mean( do.call("rbind" , lapply(iCF_D4_BT, function(df) nrow(df))) ), 0)
    mean_node_D3 <-round(mean( do.call("rbind" , lapply(iCF_D3_BT, function(df) nrow(df))) ), 0)
    mean_node_D2 <-round(mean( do.call("rbind" , lapply(iCF_D2_BT, function(df) nrow(df))) ), 0)
    
    
    #best tree list format
    iCF_D5_BT_L <- GET_TREE_L(iCF_D5, 2) 
    iCF_D4_BT_L <- GET_TREE_L(iCF_D4, 2) 
    iCF_D3_BT_L <- GET_TREE_L(iCF_D3, 2)
    iCF_D2_BT_L <- GET_TREE_L(iCF_D2, 2)
    
    
    #for split frequency, doesn't consider split values		
    iCF_D5_SF <- GET_TREE_L(iCF_D5, 3)
    iCF_D4_SF <- GET_TREE_L(iCF_D4, 3)
    iCF_D3_SF <- GET_TREE_L(iCF_D3, 3)
    iCF_D2_SF <- GET_TREE_L(iCF_D2, 3)
    
    #HTE P value for ALL iteractions of CF
    HTE_P_iCF_D5 <- HTE_P_extract(iCF_D5_BT, "D5")
    HTE_P_iCF_D4 <- HTE_P_extract(iCF_D4_BT, "D4")
    HTE_P_iCF_D3 <- HTE_P_extract(iCF_D3_BT, "D3")
    HTE_P_iCF_D2 <- HTE_P_extract(iCF_D2_BT, "D2")
    
    HTE_P_D5_median <- 	median(HTE_P_iCF_D5$HTE_P_cf)
    HTE_P_D4_median <- 	median(HTE_P_iCF_D4$HTE_P_cf)
    HTE_P_D3_median <- 	median(HTE_P_iCF_D3$HTE_P_cf)
    HTE_P_D2_median <- 	median(HTE_P_iCF_D2$HTE_P_cf)
    
  } else if (iterationNo==1){
    
    iCF_D5 <- oneCF(leafsize$D5,  treeNo, "none", "D5", split_val_round_posi) 
    iCF_D4 <- oneCF(leafsize$D4,  treeNo, "none", "D4", split_val_round_posi) 
    iCF_D3 <- oneCF(leafsize$D3,  treeNo, "none", "D3", split_val_round_posi) 
    iCF_D2 <- oneCF(leafsize$D2,  treeNo, "none", "D2", split_val_round_posi)  
    
    #tree b dataframe	
    iCF_D5_BT <- iCF_D5$besttreelist  
    iCF_D4_BT <- iCF_D4$besttreelist              
    iCF_D3_BT <- iCF_D3$besttreelist              
    iCF_D2_BT <- iCF_D2$besttreelist 
    #for split frequency
    iCF_D5_SF <- list(iCF_D5$d)  
    iCF_D4_SF <- list(iCF_D4$d)          
    iCF_D3_SF <- list(iCF_D3$d)           
    iCF_D2_SF <- list(iCF_D2$d)
  }
  #drop other column not related to tree structure
  iCF_D5_t  <-PRE_MAJORITY_TREE  (iCF_D5_BT)
  iCF_D4_t  <-PRE_MAJORITY_TREE  (iCF_D4_BT)
  iCF_D3_t  <-PRE_MAJORITY_TREE  (iCF_D3_BT)
  iCF_D2_t  <-PRE_MAJORITY_TREE  (iCF_D2_BT)
  #drop spliting value column
  iCF_D5_t_r<-PRE_MAJORITY_TREE_R(iCF_D5_BT)
  iCF_D4_t_r<-PRE_MAJORITY_TREE_R(iCF_D4_BT)
  iCF_D3_t_r<-PRE_MAJORITY_TREE_R(iCF_D3_BT)
  iCF_D2_t_r<-PRE_MAJORITY_TREE_R(iCF_D2_BT)
  #Orange column in iCF algorithm
  vote_D5_tree       <-MAJORITY_VOTE(iCF_D5_BT,         iCF_D5_t,        iCF_D5_t,        iCF_D5_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  vote_D4_tree       <-MAJORITY_VOTE(iCF_D4_BT,         iCF_D4_t,        iCF_D4_t,        iCF_D4_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  vote_D3_tree       <-MAJORITY_VOTE(iCF_D3_BT,         iCF_D3_t,        iCF_D3_t,        iCF_D3_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  vote_D2_tree       <-MAJORITY_VOTE(iCF_D2_BT,         iCF_D2_t,        iCF_D2_t,        iCF_D2_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  #Vote by key strucutre of tree (igorning splitting values)
  vote_D5_tree_R  <-MAJORITY_VOTE(iCF_D5_BT,      iCF_D5_t_r,   iCF_D5_t_r,   iCF_D5_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  vote_D4_tree_R  <-MAJORITY_VOTE(iCF_D4_BT,      iCF_D4_t_r,   iCF_D4_t_r,   iCF_D4_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  vote_D3_tree_R  <-MAJORITY_VOTE(iCF_D3_BT,      iCF_D3_t_r,   iCF_D3_t_r,   iCF_D3_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  vote_D2_tree_R  <-MAJORITY_VOTE(iCF_D2_BT,      iCF_D2_t_r,   iCF_D2_t_r,   iCF_D2_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  #obtain the synthetic (i.e. averaging splitting value among those trees with the most populuar key structure) tree structure for subgroup decision
  vote_D5_tree.syn     = vote_D5_tree_R$majority.syn
  vote_D4_tree.syn     = vote_D4_tree_R$majority.syn
  vote_D3_tree.syn.ori = vote_D3_tree_R$majority.syn
  vote_D2_tree.syn.ori = vote_D2_tree_R$majority.syn
  
  HTE_P_D5_M_median <- vote_D5_tree$HTE_P_cf_M_median
  HTE_P_D4_M_median <- vote_D4_tree$HTE_P_cf_M_median
  HTE_P_D3_M_median <- vote_D3_tree$HTE_P_cf_M_median
  HTE_P_D2_M_median <- vote_D2_tree$HTE_P_cf_M_median
  
  ###############################################
  # original split value without imputation
  ###############################################
  
  #Orange column for iCF algorithm: convert synthetic tree structure to 
  vote_D5_subgroup <- TREE2SUBGROUP(vote_D5_tree.syn )
  vote_D4_subgroup <- TREE2SUBGROUP(vote_D4_tree.syn )
  vote_D3_subgroup <- TREE2SUBGROUP(vote_D3_tree.syn.ori )
  vote_D2_subgroup <- TREE2SUBGROUP(vote_D2_tree.syn.ori )
  
  
  
  
  s_SGdeci=Sys.time()
  #White column for iCF algorihtm figure B, formular of models to predict Y* at each depth, for training and testing set, respectively.
  #has nothing to group lasso anymore, get data and models from G2-G5
  
  DATA4grplasso_iCF_tr     <- DATA4grplasso(vote_D5_subgroup, vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup, Train_ID)
  
  #according to Algorithm, we get G_D from training set, then buidl model for Y* in testing set!!!
  DATA4grplasso_iCF_te     <- DATA4grplasso(vote_D5_subgroup, vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup, Test_ID )
  
  #extracting training or testing set with g_D, and g_D definition (defined by spliting variable and splitting values!)
  Train_ID_SG_iCF       <- DATA4grplasso_iCF_tr$Dat_ID_SG
  Test_ID_SG_iCF        <- DATA4grplasso_iCF_te$Dat_ID_SG
  
  #the following formulas work for both binary and continuous outcome

  formula_basic<<-  as.formula(  paste0("Y ~ W +", paste0(colnames(X), collapse = " + ") )   )

  #main effect + interactions of W and G
  formula_g2   <<- DATA4grplasso_iCF_tr$formula_g2
  formula_g3   <<- DATA4grplasso_iCF_tr$formula_g3
  formula_g4   <<- DATA4grplasso_iCF_tr$formula_g4
  formula_g5   <<- DATA4grplasso_iCF_tr$formula_g5
  

  Deci_Final_iCF.tran.tr       <- CF_GROUP_DECISION(HTE_P_cf.raw, Train_ID_SG_iCF, "iCF", "transform", P_threshold)

  
  e_iCF=Sys.time()
  
  
  return(list ( vote_D5_tree.syn = vote_D5_tree.syn,
                vote_D4_tree.syn = vote_D4_tree.syn,
                vote_D3_tree.syn = vote_D3_tree.syn.ori,
                vote_D2_tree.syn =  vote_D2_tree.syn.ori,
                
                vote_D5_subgroup=vote_D5_subgroup,
                vote_D4_subgroup=vote_D4_subgroup,
                vote_D3_subgroup=vote_D3_subgroup,
                vote_D2_subgroup=vote_D2_subgroup,
                
                time_iCF_AIC =  difftime(e_iCF, s_iCF, units="secs"), 
                
                stability_D5_T       = unlist(vote_D5_tree$stability),        
                stability_D5_T_r     = unlist(vote_D5_tree_R$stability),    
                stability_D4_T       = unlist(vote_D4_tree$stability),        
                stability_D4_T_r     = unlist(vote_D4_tree_R$stability),      
                stability_D3_T       = unlist(vote_D3_tree$stability),        
                stability_D3_T_r     = unlist(vote_D3_tree_R$stability),       
                stability_D2_T       = unlist(vote_D2_tree$stability),   
                stability_D2_T_r     = unlist(vote_D2_tree_R$stability),
                
                Deci_Final_iCF.tran.tr  = Deci_Final_iCF.tran.tr ,
       
                Train_ID_SG_iCF      = Train_ID_SG_iCF,
                Test_ID_SG_iCF       = Test_ID_SG_iCF 
                
  )
  )  
}
