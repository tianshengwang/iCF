#' This function connects the majority vote (or any) deicison with the list of original tree in list format
#' @param Maj_Tree_Decision a majority voted tree decision df that need to be presented in a tree plot 
#' @param ori_TreeL  #a list of original tree list   
#' 
#' @return The original tree (lsit) that match the voted tree decision (df) 
#'
#' @export
#' 

GET_TREE_L <- function (iCF_D, listindex) {
  iCF_D_1 <- lapply(iCF_D,"[", listindex) 		#subsetting nested list to get besttree, extract the "listindex"th list from mother list
  if (listindex == 1 | listindex== 3 ){     # listindex=1 to extract besttreeelist, listindex=3 to extract splitfreqlist
    iCF_D_2 <- lapply(iCF_D_1, function(L) data.frame(L[[1]]) ) 
  } else {
    iCF_D_2 <- lapply(iCF_D_1, function(L) L[[1]] ) 
  }
  return(iCF_D_2)
}



#' This function for manually tuning MLS by inputting denominator (divide sample size into) and show plot of depth
#' @param demoninator a demoninator to devide the sample
#' @param depth depth descripion, e.g. "D4"   
#' 
#' @return The original tree (lsit) that match the voted tree decision (df) 
#'
#' @export
#' 
MinLeafSizeTune <- function(dat, denominator, treeNo, iterationNo, split_val_round_posi, depth ,color){
  iCF_D<- iCF_basic(denominator,  treeNo, iterationNo, dat, depth,  split_val_round_posi)
  iCF_D_BT=GET_TREE_L(iCF_D, 1)
  #turn tree into subgroup  
  iCF_D_BT_SG        <- lapply(iCF_D_BT, function(df)   TREE2SUBGROUP(df)$majority)
  #obtain depth of each leaf
  iCF_D_BT_SG_length <- lapply(iCF_D_BT_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
  #obtain depth of a tree
  iCF_D_BT_SG_Depth  <- lapply(iCF_D_BT_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
  #obtain mean depth of the causal forest
  iCF_D_BT_SG_D_mean <- round(mean( do.call("rbind" , lapply( iCF_D_BT_SG_Depth, function(df) df$Depth[1] )) ), 0) 
  
  depth_table=as.data.frame(table(do.call("rbind" , lapply( iCF_D_BT_SG_Depth, function(df) df$Depth[1] ))))%>%
              dplyr::rename(Depth=Var1, Frequency=Freq)
  
  depth_gg <- ggplot(data= depth_table, aes(x=Depth, y=Frequency)) +
                        geom_bar(stat="identity", fill= color) +
                    theme(axis.text.x = element_text(face="bold",  size=14),
                          axis.text.y = element_text( size=10)) +
                    ggtitle(paste0("Depth distribution of ",iterationNo, " best trees generated \n from ",iterationNo, " causal forests by the proposed \n minimum-leaf-size = N/", denominator, " (sample size divided by ",denominator, ") for ",depth, " iCF")) 
                  
  

  return(list(denominator= denominator,
              depth_mean=iCF_D_BT_SG_D_mean,
              depth_gg=depth_gg))
}






#' This function find right leaf size at each depth
#' @param depth depth of a forest
#' @param initial intially tried value (the denominator of N/initial, e.g. minimum-leaf-size=N/25 for D2 tree, initial=25
#' @param upper upper level if we need to get a deeper tree 
#' @param lower lower level if we need to get a shallower tree 
#' @iCF_D_BT_exp_SG_D_mean mean depth of best trees from causal forest grown by suggested initial value
#' 
#' @return the tuned leaf size 
#'
#' @export
#' 
#' 
GRID_LEAFSIZE <- function (depth, treeNo, iterationNo, initial, iCF_D_BT_exp_SG_D_mean, tune_unit){
if (   iCF_D_BT_exp_SG_D_mean ==depth){
  #Scenario 1: suggested value is good, very lucky!
  return(initial)
} else {
  #########################################
  # STEP 1: get a list of suggested value
  #########################################
  if (  iCF_D_BT_exp_SG_D_mean < depth #need a deeper tree
        #=================
        #Scenario 2: suggested value is too small
        #minimum leaf size for D5 (N/200) is too large to split,need to have a smaller (larger than initial value 200) minimum leaf size for D5, thus will count from 200 up to 800
        #=================
  ) {
    grid_MinNodeSize_D <<- seq(initial, initial*10, by=tune_unit)
    
  } else if ( iCF_D_BT_exp_SG_D_mean > depth #need a shallower tree
              #=================
              #Scenario 3: suggested value is too large
              #minimum leaf size for D5 (N/200) is too small to grow a very shallows, D5 tree,need to have a larger (smaller than initial value 200) minimum leaf size for D5, thus will count from 200 down to 100
              #=================
  ) {
    grid_MinNodeSize_D <<- seq(initial, initial/10, by=-tune_unit)
  }
  
  #########################################
  # STEP 2: try each sample size one by one
  #########################################
  for (m in 1: length(grid_MinNodeSize_D) ){
   #get iCF by proposed minimum node size
    iCF_D_experiment <- iCF_basic(grid_MinNodeSize_D[m],  treeNo, iterationNo, Ntrain, paste0("D",depth), 1)
   #get df of tree
    iCF_D_BT_experiment <- GET_TREE_L(iCF_D_experiment, 1) 
   #turn tree into subgroup  
    iCF_D_BT_experiment_SG        <- lapply(iCF_D_BT_experiment, function(df)   TREE2SUBGROUP(df)$majority)
   #obtain depth of each leaf
    iCF_D_BT_experiment_SG_length <- lapply(iCF_D_BT_experiment_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
   #obtain depth of a tree
    iCF_D_BT_experiment_SG_Depth  <- lapply(iCF_D_BT_experiment_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
   #obtain mean depth of the causal forest
    iCF_D_BT_experiment_SG_D_mean <- round(mean( do.call("rbind" , lapply( iCF_D_BT_experiment_SG_Depth, function(df) df$Depth[1] )) ), 0) 
    
    
    if ( iCF_D_BT_experiment_SG_D_mean == depth  )     break #if mean node # of best trees from iCF equal to 3, i.e. D2 tree structure
    D_tune = grid_MinNodeSize_D[m]
    
  }
  
  return(D_tune)
  
} 
}
  

#' This function tune leaf size for iCF 
#' @param Ntrain sample size
#' @param treeNo  tree No, need to be >= 50, if =10 has this error: Error in lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,  : 0 (non-NA) cases
#' @param iterationNo  iteration No  
#' 
#' @return the tuned leaf size 
#'
#' @export
#' 
LEAFSIZE_tune <- function(Ntrain, initial_D2, treeNo, iterationNo, tune_unit){
  
  #suggested N to start with, then check, if doesn't work then tune by increasing or  minimum leafsize
  leafsize_initial = list(D5=initial_D2*4,
                          D4=initial_D2*3, 
                          D3=initial_D2*2,
                          D2=initial_D2)
  
  #=================
  #CHECK: 
  #Scenario 1: suggested value works! if node # of D2 tree==3, then return leafsize, otherwise rerun tuning!
  #=================
  #########################
  # Checking
  #########################
  iCF_D2_exp = iCF_basic(leafsize_initial$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi)
  iCF_D2_BT_exp <<- GET_TREE_L(iCF_D2_exp, 1) 
  iCF_D2_BT_exp_D <- lapply(iCF_D2_BT_exp, function(df) nrow(df)   )
  
  #the initial_D2 (min.node.size) for D=2 tree may be too large so that the tree don't split at all, 
  #i.e. there is only ONE node-01, node-02 and node-03 don't exist at all.
  if (  min( unlist(iCF_D2_BT_exp_D) ) == 1 ) {

    leafsize_initial = list(D5=round(initial_D2*7,0), 
                            D4=round(initial_D2*5,0),  
                            D3=round(initial_D2*4,0), 
                            D2=round(initial_D2*2,0)
                            )
    
    iCF_D2_exp = iCF_basic(leafsize_initial$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi)
    iCF_D2_BT_exp <<- GET_TREE_L(iCF_D2_exp, 1) 
  }

  
  
  iCF_D3_exp = iCF_basic(leafsize_initial$D3,  treeNo, iterationNo, Ntrain, "D3",  split_val_round_posi)
  iCF_D3_BT_exp <<- GET_TREE_L(iCF_D3_exp, 1) 
  
  iCF_D4_exp = iCF_basic(leafsize_initial$D4,  treeNo, iterationNo, Ntrain, "D4",  split_val_round_posi)
  iCF_D4_BT_exp <<- GET_TREE_L(iCF_D4_exp, 1) 

  iCF_D5_exp = iCF_basic(leafsize_initial$D5,  treeNo, iterationNo, Ntrain, "D5",  split_val_round_posi)
  iCF_D5_BT_exp <<- GET_TREE_L(iCF_D5_exp, 1) 
  
  iCF_D2_BT_exp_SG        <- lapply(iCF_D2_BT_exp, function(df)   TREE2SUBGROUP(df)$majority)
  iCF_D2_BT_exp_SG_length <- lapply(iCF_D2_BT_exp_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
  iCF_D2_BT_exp_SG_Depth  <- lapply(iCF_D2_BT_exp_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
  iCF_D2_BT_exp_SG_D_mean <- round(mean( do.call("rbind" , lapply( iCF_D2_BT_exp_SG_Depth, function(df) df$Depth[1] )) ), 0) 
  
  
  iCF_D3_BT_exp_SG        <- lapply(iCF_D3_BT_exp, function(df)   TREE2SUBGROUP(df)$majority)
  iCF_D3_BT_exp_SG_length <- lapply(iCF_D3_BT_exp_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
  iCF_D3_BT_exp_SG_Depth  <- lapply(iCF_D3_BT_exp_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
  iCF_D3_BT_exp_SG_D_mean <- round(mean( do.call("rbind" , lapply( iCF_D3_BT_exp_SG_Depth, function(df) df$Depth[1] )) ), 0) 
  
  
  iCF_D4_BT_exp_SG        <- lapply(iCF_D4_BT_exp, function(df)   TREE2SUBGROUP(df)$majority)
  iCF_D4_BT_exp_SG_length <- lapply(iCF_D4_BT_exp_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
  iCF_D4_BT_exp_SG_Depth  <- lapply(iCF_D4_BT_exp_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
  iCF_D4_BT_exp_SG_D_mean <- round(mean( do.call("rbind" , lapply( iCF_D4_BT_exp_SG_Depth, function(df) df$Depth[1] )) ), 0) 
  
  
  iCF_D5_BT_exp_SG        <- lapply(iCF_D5_BT_exp, function(df)   TREE2SUBGROUP(df)$majority)
  iCF_D5_BT_exp_SG_length <- lapply(iCF_D5_BT_exp_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
  iCF_D5_BT_exp_SG_Depth  <- lapply(iCF_D5_BT_exp_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
  iCF_D5_BT_exp_SG_D_mean <- round(mean( do.call("rbind" , lapply( iCF_D5_BT_exp_SG_Depth, function(df) df$Depth[1] )) ), 0) 
  
  
 # D5_tune <- GRID_LEAFSIZE(5, treeNo,  iterationNo, leafsize_initial$D5, iCF_D5_BT_exp_SG_D_mean, tune_unit)
  
  if ( is.null( tryCatch( GRID_LEAFSIZE(5, treeNo,  iterationNo, leafsize_initial$D5, iCF_D5_BT_exp_SG_D_mean, tune_unit), error=function(e){}) ) == TRUE ) {
    warning("Failed to find D5_tune value, use 175")
    D5_tune = 175
  } else {
    D5_tune <-            GRID_LEAFSIZE(5, treeNo,  iterationNo, leafsize_initial$D5, iCF_D5_BT_exp_SG_D_mean, tune_unit)
  }
#  D4_tune <- GRID_LEAFSIZE(4, treeNo,  iterationNo, leafsize_initial$D4, iCF_D4_BT_exp_SG_D_mean, tune_unit)
  if ( is.null( tryCatch( GRID_LEAFSIZE(4, treeNo,  iterationNo, leafsize_initial$D4, iCF_D4_BT_exp_SG_D_mean, tune_unit), error=function(e){}) ) == TRUE ) {
    warning("Failed to find D4_tune value, use 125")
    D4_tune = 125 
  } else {
    D4_tune <-            GRID_LEAFSIZE(4, treeNo,  iterationNo, leafsize_initial$D4, iCF_D4_BT_exp_SG_D_mean, tune_unit)
  }
  
  #D3_tune <- GRID_LEAFSIZE(3, treeNo,  iterationNo, leafsize_initial$D3, iCF_D3_BT_exp_SG_D_mean, tune_unit)
  if ( is.null( tryCatch( GRID_LEAFSIZE(3, treeNo,  iterationNo, leafsize_initial$D3, iCF_D3_BT_exp_SG_D_mean, tune_unit), error=function(e){}) ) == TRUE ) {
    warning("Failed to find D3_tune value, use 100")
    D3_tune = 100 
  } else {
    D3_tune <-            GRID_LEAFSIZE(3, treeNo,  iterationNo, leafsize_initial$D3, iCF_D3_BT_exp_SG_D_mean, tune_unit)
  }
  
  
  #D2_tune <- GRID_LEAFSIZE(2, treeNo,  iterationNo, leafsize_initial$D2, iCF_D2_BT_exp_SG_D_mean, tune_unit)
  if ( is.null( tryCatch( GRID_LEAFSIZE(2, treeNo,  iterationNo, leafsize_initial$D2, iCF_D2_BT_exp_SG_D_mean, tune_unit), error=function(e){}) ) == TRUE ) {
    warning("Failed to find D2_tune value, use 50")
    D2_tune = 50 
  } else {
    D2_tune <-            GRID_LEAFSIZE(2, treeNo,  iterationNo, leafsize_initial$D2, iCF_D2_BT_exp_SG_D_mean, tune_unit)
  }
  
    #########################
    # STEP 3. double check
    #########################
    leafsize_tune <- list(D5=D5_tune, #takes no computation time
                          D4=D4_tune, #takes no computation time
                          D3=D3_tune, 
                          D2=D2_tune)
   return(leafsize_tune)
}






#' This function is the wrapper function to run raw causal forest
#' @param Train_cf training data
#' @param min.split.var minimum number of variables to split population
#' 
#' @return the key results from raw causal forest 
#'
#' @export
#' 

CF_RAW_key <- function(Train_cf, min.split.var, variable_type, hdpct){
  s_rawCF =  Sys.time()
  #------------------------------------------------------------------------
  #step 1. get W.hat, Y.hat based on X
  #------------------------------------------------------------------------
  X <- Train_cf[,vars_forest]
  Y <- as.vector( as.numeric( Train_cf[,1] ) )
  W <- as.vector( as.numeric( Train_cf[,2] ) )  

  W.forest <- regression_forest(X, W)
  W.hat    <- predict(W.forest)$predictions
  
  #forest-based method to predict Y
  Y.forest <- regression_forest(X, Y)
  Y.hat    <- predict(Y.forest)$predictions
  
  #------------------------------------------------------------------------
  #step 2. run CF to get VI and P-value of omnibus test for HTE presence 
  #------------------------------------------------------------------------
  cf.raw = causal_forest(X, Y,  W,  Y.hat = Y.hat, W.hat = W.hat)
  
  best_tree_info<-find_best_tree(cf.raw, "causal")
  best_tree_info$best_tree
  # Plot trees
  par(mar=c(1,1,1,1))
  tree.plot = plot(grf::get_tree(cf.raw, best_tree_info$best_tree))
  tree.plot
  
  HTE_P_cf.raw <-test_calibration(cf.raw)[8]
  HTE_P_cf.raw
  varimp_cf <- variable_importance(cf.raw)
  #par(ask=F)
  #old.par <- par(mar = c(0, 0, 0, 0))
  #par(old.par)
  #par()
  plot(varimp_cf)
  text(1:ncol(X), varimp_cf,labels=colnames(X))
  
  #---------------------------------------------------------------
  #if usng > mean critiera, may select only 1 varialbe when using "binary" (indicator) ,
  selected_cf.idx_mean = which(varimp_cf > mean(varimp_cf)) 
  selected_cf.idx_q3 = which(varimp_cf > quantile(varimp_cf, 0.75) )
  selected_cf.idx_q2 = which(varimp_cf > quantile(varimp_cf, 0.5) )
  selected_cf.idx_q1 = which(varimp_cf > quantile(varimp_cf, 0.25) )
  selected_cf.idx_q0 = which(varimp_cf > 0 ) #selected all
  selected_cf.idx_hd = which(varimp_cf > quantile(varimp_cf, hdpct) ) #for hdiCF
  #to develop D5 tree, need to select >=4 variables. e.g., it will develop D4 tree if selecting only 3 variables.
  #28=7*4, i.e., 7 is the minimum node# for a depth 4symmetrical tree, and when 7 variables are the q4 part of distribution, it requires 28 varaibles  
  #60=15*4, i.e., 15 is the minimum node# for a depth 5symmetrical tree, and when 7 variables are the q4 part of distribution, it requires 28 varaibles 
  if(variable_type=="hd"){
  selected_cf.idx =  selected_cf.idx_hd 
  } else if (variable_type=="hdshrink") {
  selected_cf.idx =  selected_cf.idx_q0   
  } else{
  if( length(selected_cf.idx_mean) >=  min.split.var ) {
    selected_cf.idx = selected_cf.idx_mean
  } else if(  length(selected_cf.idx_mean) <  min.split.var &  length( selected_cf.idx_q3) >=  min.split.var
  ) {
    selected_cf.idx = selected_cf.idx_q3
  } else if(  length(selected_cf.idx_mean) <  min.split.var & length( selected_cf.idx_q3) <  min.split.var &  length( selected_cf.idx_q2) >=  min.split.var
  ) {
    selected_cf.idx = selected_cf.idx_q2
  } else if(  length(selected_cf.idx_mean) <  min.split.var & length( selected_cf.idx_q2) <  min.split.var & length( selected_cf.idx_q1) >=  min.split.var 
  ) {
    selected_cf.idx = selected_cf.idx_q1
  } else if (        length(selected_cf.idx_mean) <  min.split.var & length( selected_cf.idx_q1) < min.split.var #& 
     # dim(X)[2]>=15*4 # number of covariates >= 60
  ) {
    selected_cf.idx =  selected_cf.idx_q0 # select all
  } 
  }
  #---------------------------------------------------------------
  #head( X[, selected_cf.idx], 1 ) #be careful!!! not head(Train[, selected_cf.idx])
  
  #colnames(X[, selected_cf.idx ] )
  e_rawCF =  Sys.time()
  time_rawCF =  difftime(e_rawCF, s_rawCF,units="secs")
  
  return(list(X=X, 
              Y=Y,
              W=W,
              Y.hat=Y.hat,
              W.hat=W.hat,
              HTE_P_cf.raw=HTE_P_cf.raw,
              varimp_cf=varimp_cf,
              selected_cf.idx=selected_cf.idx,
              time_rawCF=time_rawCF))
  
}

#' This function is the wrapper function to run raw causal forest and for initial variable selection for hdICF, 
#' which will be applied to a list of selected X by incremetally trimming VI score 
#' @param Train_cf training data
#' @param X_selected selected variables index
#' 
#' @return the accuracy of predicting W from raw causal forest with different selected hdX
#'
#' @export
#' 
#Xname_selected[[19]]

hdiCF_VS_incremental <- function(#Train_cf, 
                                 X_selected){
  #X <- Train_cf[, X_selected]
  X <- X[,X_selected]
  
  #Y <- as.vector( as.numeric( Train_cf[,1] ) )
  #W <- as.vector( as.numeric( Train_cf[,2] ) )  
  
  #X_untreat <- Train_cf[, X_selected][ which(W==0),]
  #Y_untreat <- as.vector( as.numeric( Train_cf[ which(W==0),][,1] ) )
  #W_untreat <- as.vector( as.numeric( Train_cf[ which(W==0),][,2] ) )  

  #------------------------------------------------------------------------
  #when incrementally trim tail of VI score from all varaibles (VI doesn't change) only need W.hat 
  W.forest <- regression_forest(X, W)
  W.hat    <- predict(W.forest)$predictions
  
  #forest-based method to predict Y
 # Y.forest_untreat <- regression_forest(X_untreat, Y_untreat)
 # Y.hat_untreat    <- predict(Y.forest_untreat)$predictions
  
  # Compute confusion matrix
  
 CMatrix_W <- caret::confusionMatrix(data = as.factor(ifelse(W.hat>0.5,1,0)) ,
                                  reference = as.factor(W),
                                  positive = "1")
 accuracy_W <-  as.numeric(CMatrix_W$overall[1])
 
 #doesn't make much sense using only X (without W) to predict Y
#CMatrix_Y_untreat <- caret::confusionMatrix(data = as.factor(ifelse(Y.hat_untreat>0.5,1,0)) ,
#                                     reference = as.factor(Y_untreat),
#                                    positive = "1")
#accuracy_Y_untreat <-  as.numeric(CMatrix_Y_untreat$overall[1])
 
  #create transformed outcome
 # Y_star <- ifelse(W==1, Y/W.hat, -Y/(1-W.hat))
  
  #forest-based method to predict Y_star
  #Y.forest.star <- regression_forest(X, Y_star)
  #Y.hat.star    <- predict(Y.forest.star)$predictions

  
  #RMSE_Y_star = sqrt( sum((Y_star - Y.hat.star)^2)/length(Y)) 
return(#list(
       accuracy_W = accuracy_W#, 
       #accuracy_Y_untreat = accuracy_Y_untreat)
       )
}
  
 


hdiCF_VS_iterative <- function(#Train_ori, 
  s, cutoffpct){
  accuracy_W_list = list()
  X_ncol_list = list()
  CMatrix_W_list = list()
  HTE_P_cf_list = list()
  Train_list = list()  
  Xselectedname_list = list()
  
  for (s in 1:s) { 
    if(s==1){
      #All HD variables, original dataset 
      # Train <<- Train_ori
      W.forest <- regression_forest(X, W)
      W.hat    <- predict(W.forest)$predictions
      
      Y.forest <- regression_forest(X, Y)
      Y.hat    <- predict(Y.forest)$predictions
      
      cf.raw = causal_forest(X, Y,  W,  Y.hat = Y.hat, W.hat = W.hat)
      HTE_P_cf.raw <-test_calibration(cf.raw)[8]
      varimp_cf <- variable_importance(cf.raw)
      selected_cf.idx <<- which(varimp_cf >= quantile(varimp_cf, cutoffpct ) )
      
    } else {
      #iteratively generate new dataset, new raw CF from new W.hat and new Y.hat, new VI score
      X <<- X[,selected_cf.idx] #from previous (s-1)th loop
      
      W.forest <- regression_forest(X, W)
      W.hat    <- predict(W.forest)$predictions
      #new Y.hat
      Y.forest <- regression_forest(X, Y)
      Y.hat    <- predict(Y.forest)$predictions
      
      cf.raw = causal_forest(X,   Y,  W,  Y.hat = Y.hat, W.hat = W.hat)
      HTE_P_cf.raw <-test_calibration(cf.raw)[8]
      varimp_cf <- variable_importance(cf.raw) 
      selected_cf.idx <<- which(varimp_cf >= quantile(varimp_cf, cutoffpct ) ) #for next (s+1)th loop, update it
      
    }
    #X <- Train_cf[, X_selected]
    #cf_raw_key.tr <- CF_RAW_key(Train, 1, "hd", hdpct=cutoffpct) 
    # varimp_cf  <- cf_raw_key.tr$varimp_cf
    #  selected_cf.idx = cf_raw_key.tr$selected_cf.idx #which(varimp_cf >= quantile(varimp_cf, cutoffpct) )
    #  length(selected_cf.idx)
    #W.hat <-   cf_raw_key.tr$W.hat
    Xselectedname_list [[s]] <- colnames(X[,c(selected_cf.idx)])
    
    CMatrix_W <- caret::confusionMatrix(data = as.factor(ifelse(W.hat>0.5,1,0)) ,
                                        reference = as.factor(W),
                                        positive = "1")
    
    CMatrix_W_list[[s]] <- CMatrix_W
    HTE_P_cf_list[[s]]= HTE_P_cf.raw
    accuracy_W_list[[s]] <-  as.numeric(CMatrix_W$overall[1])
    
    #redefine X, Train
    # X <<- X[,c(selected_cf.idx)] #<<- to make it avaialbe outside of a specific loop
    # ncol(X)
    #  Train <<- Train[,c("Y", "W", colnames( X )) ]
    #vars_forest <<- colnames( Train %>% dplyr::select(-c("Y", "W")))
    X_ncol_list[[s]] <- ncol(X)
    
    #Train_list[[s]] <- Train
    
  }
  all_list <- rlist::list.zip(accuracy_W_list, CMatrix_W_list,  X_ncol_list, Xselectedname_list, HTE_P_cf_list #, Train_list
  )
  return(all_list) 
  
}




CONVERT_tree_IT2CF <- function(btree_it) {
  btree_cf <- btree_it  %>%  
    mutate (is_leaf        = ifelse(is.na(vname)==TRUE , TRUE , FALSE))  %>%
    mutate (left_child     = ifelse(is_leaf=="FALSE", paste0(node,1), "NA")) %>%
    mutate (right_child    = ifelse(is_leaf=="FALSE", paste0(node,2), "NA")) %>%
    mutate (split_variable = ifelse(is.na(vname)==FALSE, paste0(as.character(vname)), "NA" ) ) %>%
    mutate (split_value    = ifelse(is.na(best.cut)==FALSE, round(as.numeric(as.character(best.cut)),1), "NA") ) %>%
    mutate (samples        = ifelse(is_leaf=="TRUE", n, "NA") ) %>%
    mutate (avg_Y          = 0) %>%
    mutate (avg_W          = 0) %>%
    mutate (b              = 1) %>%
    mutate (HTE_P_cf       = 0) %>%
    subset (select=  c(node, is_leaf, left_child, right_child, split_variable, split_value, samples, avg_Y, avg_W, b, HTE_P_cf)) %>%
    mutate (node = paste("node", node, sep = '-') ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -2, -2) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -3, -3) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -4, -4) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -5, -5) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (left_child  = ifelse(nchar(left_child)  == 2 & left_child  !="NA", paste0("000", left_child),  left_child ) ) %>%
    mutate (left_child  = ifelse(nchar(left_child)  == 3 & left_child  !="NA", paste0("00" , left_child),  left_child ) ) %>%
    mutate (left_child  = ifelse(nchar(left_child)  == 4 & left_child  !="NA", paste0("0" ,  left_child),  left_child ) ) %>%
    mutate (right_child = ifelse(nchar(right_child) == 2 & right_child !="NA", paste0("000", right_child), right_child ) ) %>%
    mutate (right_child = ifelse(nchar(right_child) == 3 & right_child !="NA", paste0("00" , right_child), right_child ) ) %>%
    mutate (right_child = ifelse(nchar(right_child) == 4 & right_child !="NA", paste0("0" ,  right_child), right_child ) ) %>%
    arrange(node) 
  return( as.data.frame(btree_cf) )       
}


CONVERT_tree_stima2CF <- function(btree_stima) {
btree_cf <- btree_stima %>% 
mutate (node     = 1:nrow(btree_stima)) %>%
mutate (node     = paste("node",node, sep = "-" )) %>%
mutate (node = ifelse( stringr::str_sub(node,  -2, -2) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
mutate (is_leaf     = ifelse(nchar(Region)!=0 , TRUE , FALSE)) %>%
mutate ( value_1    =  lead(Splitpoint,n=1)   )%>%  
mutate ( value_2    =  lead(Splitpoint,n=2)   )%>% 
mutate ( value_3    =  lead(Splitpoint,n=3)   )%>%  
mutate ( value_4    =  lead(Splitpoint,n=4)   )%>% 
mutate ( var_1      =  lead(Predictor,n=1)   )   %>%  
mutate ( var_2      =  lead(Predictor,n=2)   )   %>%  
mutate ( var_3      =  lead(Predictor,n=3)   )   %>%  
mutate ( var_4      =  lead(Predictor,n=4)   )   %>% 
mutate ( Sign_1          =  lead(Sign,n=1)   )   %>%  
mutate ( Sign_2          =  lead(Sign,n=2)   )   %>%  
mutate ( Sign_3          =  lead(Sign,n=3)   )   %>%  
mutate ( Sign_4          =  lead(Sign,n=4)   )   %>% 
mutate ( node_1     =  lead(node,n=1)   )   %>%  
mutate ( node_2     =  lead(node,n=2)   )   %>%  
mutate ( node_3     =  lead(node,n=3)   )   %>%  
mutate ( node_4     =  lead(node,n=4)   )   %>% 
mutate ( Region_1     =  lead(Region,n=1)   )   %>%  
mutate ( Region_2     =  lead(Region,n=2)   )   %>%  
mutate ( Region_3     =  lead(Region,n=3)   )   %>%  
mutate ( Region_4     =  lead(Region,n=4)   )   %>% 
mutate (left_child =
ifelse(is_leaf=="TRUE", "NA",          
ifelse((Sign_1 == "<=" & value_1== value_2 & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, stringr::str_sub(node_1,-2,-1) , 
ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                        & value_2== value_3 & var_2==var_3 &  ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, stringr::str_sub(node_2,-2,-1),
ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                        & value_1== value_2 & #value_3== value_4 & 
                                                            var_1==var_2 #& var_3==var_4 
                                                            & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, stringr::str_sub(node_3,-2,-1), NA
))))) %>%
mutate (right_child =
ifelse(is_leaf=="TRUE", "NA",
ifelse((Sign_1 == "<=" & value_1== value_2 & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, stringr::str_sub(node_2,-2,-1) , 
ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                        & value_2== value_3 & var_2==var_3 & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) ) ==TRUE, stringr::str_sub(node_3,-2,-1),
ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                        & value_1== value_2 & value_3== value_4 & var_1==var_2 & var_3==var_4 & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, stringr::str_sub(node_4,-2,-1), NA
))))) %>%
mutate (split_variable = 
ifelse(is_leaf=="TRUE", "NA",          
ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"  
                                        & value_1== value_2 #& value_3== value_4 
                                                            &( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, var_1, 
ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                        #& value_2== value_3 & var_2==var_3 
                                        & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) ) ==TRUE,    var_2,
ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                        & value_1== value_2 & value_3== value_4 & var_1==var_2 & var_3==var_4 & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, var_3, NA
))))) %>%
mutate (split_value = 
ifelse(is_leaf=="TRUE", "NA",          
ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"  
                                        & value_1== value_2 #& value_3== value_4 
                                                            & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, round(as.numeric(as.character(value_1)),1), 
ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                        #& value_2== value_3 & var_2==var_3 
                                        & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) ) ==TRUE,   round(as.numeric(as.character(value_2)),1),
ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                        & value_1== value_2 & value_3== value_4 & var_1==var_2 & var_3==var_4 & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, round(as.numeric(as.character(value_3)),1), NA
)))))  %>%

    mutate (samples        = ifelse(is_leaf=="TRUE", n, "NA") ) %>%
    mutate (avg_Y          = 0) %>%
    mutate (avg_W          = 0) %>%
    mutate (b              = 1) %>%
    mutate (HTE_P_cf       = 0) %>%
    subset (select=  c(node, is_leaf, left_child, right_child, split_variable, split_value, samples, avg_Y, avg_W, b, HTE_P_cf))

  return( as.data.frame(btree_cf) )       
}


GET_TREE_DF <- function(No_nodes, treetype, split_val_round_posi) {
  node = rep(NA, No_nodes)
  is_leaf = rep(NA,  No_nodes)
  left_child = rep(NA,  No_nodes)
  right_child = rep(NA,  No_nodes)
  split_variable = rep(NA, No_nodes)
  split_value = rep(NA, No_nodes)
  samples = rep(NA,  No_nodes)
  avg_Y = rep(NA, No_nodes)
  avg_W = rep(NA,  No_nodes)
  
  for (i in 1 : No_nodes) {
    node[i] <- paste('node', i, sep = '-')
    
    #add leading 0 to make node-01 format for number less than 10
    if(stringr::str_sub(node[i],  -2, -2) =="-" ){node[i] <- sub("-([0-9])", "-0\\1", node[i])}
    
    #if $is_leaf = FALSE, i.e. if not leaf
    if (treetype$nodes[[i]][1]=="FALSE") {
      is_leaf[i] <- "FALSE"
      left_child[i] <- treetype$nodes[[i]]$left_child
      right_child[i] <- treetype$nodes[[i]]$right_child
      
      #add leading 0 to make node-01 format for number less than 10
      if(nchar(left_child[i])  ==1 ){left_child[i]  <- sub("([0-9])", "0\\1",  left_child[i])}  else {left_child[i] =left_child[i]}
      if(nchar(right_child[i]) ==1 ){right_child[i] <- sub("([0-9])", "0\\1",  right_child[i])} else {right_child[i]=right_child[i]}
      
      
      split_variable[i] <- treetype$columns[treetype$nodes[[i]]$split_variable]  
      split_value[i] <- round( treetype$nodes[[i]]$split_value, split_val_round_posi)                
      samples[i] <-"NA"
      avg_Y[i] <-"NA"
      avg_W[i] <-"NA"
    } else {
      #if $is_leaf = TRUE, i.e. if leaf 
      is_leaf[i] <- "TRUE"
      left_child[i] <- "NA"
      right_child[i] <- "NA"
      split_variable[i] <- "NA"  
      split_value[i] <- "NA"
      samples[i] <-length(treetype$nodes[[i]]$samples)
      avg_Y[i] <-treetype$nodes[[i]]$leaf_stats[1]
      avg_W[i] <-treetype$nodes[[i]]$leaf_stats[2]
    }
  }
  
  tree_information <- data.frame(node, is_leaf, left_child, right_child, split_variable, split_value, samples, avg_Y, avg_W, stringsAsFactors = F)
  return(tree_information)
}

PLOT_BT<-function(forest, type){
best_tree_info<-find_best_tree(forest, type)
best_tree_info$best_tree
plot_bt<-plot(grf::get_tree(forest, best_tree_info$best_tree))
return(list(best_tree_info$best_tree, plot_bt))
}




#' This function connects the majority vote (or any) deicison with the list of original tree in list format
#' @param Maj_Tree_Decision a majority voted tree decision df that need to be presented in a tree plot 
#' @param ori_TreeL  #a list of original tree list   
#' 
#' @return The original tree (lsit) that match the voted tree decision (df) 
#'
#' @export
#' 
#' 
GET_ORI_TREE_L <- function(Maj_Tree_Decision,ori_TreeL){
 Tree_PreMaj <- PRE_MAJORITY_TREE(ori_TreeL)   #a list of tree df that remove samples, avg_Y, avg_W, k/b, HTE_P_cf #k=iteraction #
  if_majority <-lapply(Tree_PreMaj , function(df) identical(df, Maj_Tree_Decision))  #identical or not label
#STANDARD_CHECK: This function connects if identical or not label to a list and subsetting the list to exclude those not identical to reference list
tree_ID_list <- STANDARD_CHECK(if_majority,  #identical or not label
                               ori_TreeL #a list of original tree list
                              )
#randomly pick a tree with the same structure as voted tree decision, as "samples, avg_Y, avg_W" doesn't mattter, these columns varied across trees due to random subsampling
tree_ID_random <- sample(tree_ID_list, 1)[[1]] 
return(tree_ID_random)
}


HTE_P_extract <- function(iCF_D, tree_depth){
  HTE_P_cf_list <- lapply(iCF_D, function(df) mean(df$HTE_P_cf))
  HTE_P_cf   <- as.data.frame(rlist::list.rbind(HTE_P_cf_list)) #Bind all list elements by row
  colnames(HTE_P_cf) <- c("HTE_P_cf")
  HTE_P_cf$tree_depth <- tree_depth
  HTE_P_cf <- tibble::rowid_to_column(HTE_P_cf, "k")
  return(HTE_P_cf)
}

HTE_P_extract_M <- function(iCF_D_M, tree_depth){
  HTE_P_cf_M <- iCF_D_M$HTE_P_cf_M
  HTE_P_cf_M$tree_depth <- tree_depth
  return(HTE_P_cf_M)
}

SF_extract_M <- function(iCF_D_M, tree_depth){
  SF_M <- iCF_D_M$SF_ID_LIST
  return(SF_M)
}

#combine all df function:
COMBO_SAVE <- function(List,  Listname){
  allonelist = do.call(rbind, List)
  #record allbesttree data
  write.table(allonelist,file=paste(Listname, depth, treeNo, iterationNo, Ntrain, tree_depth, intTRUE, ".csv", sep="_"),row.names=FALSE,col.names=T,sep=",",append=TRUE)
  #save the list of frames
  save(List,             file=paste(Listname, depth, treeNo, iterationNo, Ntrain, tree_depth, intTRUE, ".Rda", sep="_"))
}

FOREST_LIST <- function(dataset) {
  vars_forest =c("X1","X2","X3","X4","X5","X6","X8","X9","X10") 
  vars_IV=c("X7")
  X<-dataset[,vars_forest]
  Y<-dataset[,"Y"]
  W<-dataset[,"a"]
  Z<-dataset[,vars_IV]  
  forest.list <- list(X = X, Y = Y, W = W, Z = Z)
  return(forest.list)	    
}


#' This function find out a splitter is continuous, ordinal, orbinary variable
#' @param tree the dataframe for a causal tree
#' @param row_index  #the row number of the dataframe   
#' 
#' @return The level of this splitter, e.g.if continous splitter, then > 8 
#'
#' @export
#' 
#' 
SPLITTER_LEVEL <- function(tree, row_index){
  if (tree[row_index,"split_variable"] == "NA") {#if leaf
    level_splitter = 0
  } else {#if splitter
    level_splitter = length(unique (  eval(parse(text=paste("X$", tree[row_index,"split_variable"], sep = ""))) ))
  }
  return(level_splitter)
}

#' This function impute Impute split values for splitters of D3 voted tree from D4 synthetic tree if top 3 nodes are identical
#' @param level_D3_splitter_1 the node 1 (splitter1) of D3 voted tree structure
#' @param level_D3_splitter_2 the node 2 (splitter2) of D3 voted tree structure
#' @param level_D3_splitter_3 the node 3 (splitter3) of D3 voted tree structure
#' @param vote_D3_tree.syn.ori  the original vote D3 tree stucutre (without imputing split value)   
#' @param vote_D4_tree.syn  the original vote D4 tree stucutre (without imputing split value)   
#' @return The imputed D3 tree if Node 1,2,3 of D3 identical to Node 1,2, 3 of D4, or original D3 tree if Node 1, 2, 3 of D3 not identical to D4
#'
#' @export
#' 
#'
SPLIT_VAL_IMPUTE_D3 <- function(level_D3_splitter_1, level_D3_splitter_2, level_D3_splitter_3, vote_D3_tree.syn.ori, vote_D4_tree.syn     ){
  if ( identical(vote_D3_tree.syn.ori[1:3,1:5], 
                 vote_D4_tree.syn[1:3,1:5]) & #D3 and D4 trees have the same node 123
       max(level_D3_splitter_1, level_D3_splitter_2, level_D3_splitter_3) > 2 #one of the splitter is not binary
  ){
    #node 1, 2, 3, if continuous splitter, then impute splitting value from D4 tree:
    if (level_D3_splitter_1 > 2){vote_D3_tree.syn.ori[,"split_value"][ which( vote_D3_tree.syn.ori$node == "node-01" )]  <-  vote_D4_tree.syn[1,"split_value"] }
    if (level_D3_splitter_2 > 2){vote_D3_tree.syn.ori[,"split_value"][ which( vote_D3_tree.syn.ori$node == "node-02" )]  <-  vote_D4_tree.syn[2,"split_value"] }
    if (level_D3_splitter_3 > 2){vote_D3_tree.syn.ori[,"split_value"][ which( vote_D3_tree.syn.ori$node == "node-03" )]  <-  vote_D4_tree.syn[3,"split_value"] }
  }
  return(vote_D3_tree.syn.ori)
}

#' This function impute Impute split values for splitters of D2 voted tree from D4 synthetic tree if top 1 node, i.e.,root node, are identical
#' @param level_D2_splitter_1 the node 1 (splitter1) of D3 voted tree structure
#' @param vote_D2_tree.syn.ori  the original vote D3 tree stucutre (without imputing split value)   
#' @param vote_D4_tree.syn  the original vote D4 tree stucutre (without imputing split value)   
#' @return The imputed D3 tree if Node 1,2,3 of D3 identical to Node 1,2, 3 of D4, or original D3 tree if Node 1, 2, 3 of D3 not identical to D4
#'
#' @export
#' 
#'
SPLIT_VAL_IMPUTE_D2 <- function(level_D2_splitter_1, vote_D2_tree.syn.ori, vote_D4_tree.syn){
if ( identical(vote_D2_tree.syn.ori[1,1:5], 
               vote_D4_tree.syn[1,1:5]) & #D2 and D4 trees have the same node 1
     level_D2_splitter_1 > 2 #node 1 is not binary splitter
){
  vote_D2_tree.syn.ori[,"split_value"][ which( vote_D2_tree.syn.ori$node == "node-01" )]  <-  vote_D4_tree.syn[1,"split_value"] 
}
return(vote_D2_tree.syn.ori)  
}