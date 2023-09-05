#' GET_TREE_L 
#' 
#' This function obtains dataframe of best tree from output of the iCF function
#' @param iCF_D output of the iCF function 
#' @param listindex index to help extract results from iCF_D
#' 
#' @return a list of dataframes of the best trees
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



#' MinLeafSizeTune
#' 
#' This function tests the proposed denominator for minium leaf size (N/denominator) and provides the distribution of depth of generated best trees 
#' @param dat dataset
#' @param demoninator a demoninator to devide the sample size N
#' @param treeNo tree No
#' @param iterationNo iteration No (if = 1 then oneCF)
#' @param split_val_round_posi split value rounding position
#' @param depth depth descripion, e.g. "D4"
#' @param color the color for histogram plot for distribution of depth of generated best trees   
#' 
#' @return The mean and histogram of depth of generated best trees
#'
#' @export
#'
 
MinLeafSizeTune <- function(dat, denominator, treeNo, iterationNo, split_val_round_posi, depth ,color){
  iCF_D<- iCF(denominator,  treeNo, tunepara = "none", iterationNo, dat, depth,  split_val_round_posi)
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
  
  depth_gg <- ggplot2::ggplot(data= depth_table, aes(x=Depth, y=Frequency)) +
                        geom_bar(stat="identity", fill= color) +
                    theme(axis.text.x = element_text(face="bold",  size=14),
                          axis.text.y = element_text( size=10)) +
                    ggtitle(paste0("Depth distribution of ",iterationNo, " best trees generated \n from ",iterationNo, " causal forests by the proposed \n minimum-leaf-size = N/", denominator, " (sample size divided by ",denominator, ") for ",depth, " iCF")) 
                  
  return(list(denominator= denominator,
              depth_mean=iCF_D_BT_SG_D_mean,
              depth_gg=depth_gg))
}





#'GRID_LEAFSIZE
#'
#' This function find right leaf size at each depth
#' @param depth depth of a forest
#' @param treeNo tree No
#' @param iterationNo iteration No (if = 1 then oneCF)
#' @param initial intially tried value (the denominator of N/initial, e.g. minimum-leaf-size=N/25 for D2 tree, initial=25
#' @param iCF_D_BT_exp_SG_D_mean mean depth of best trees from causal forest grown by suggested initial value
#' @param tune_unit unit for tuning
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
  

#'LEAFSIZE_tune
#'
#' This function tunes leaf size for iCF 
#' @param Ntrain sample size
#' @param initial_D2 initial dominator used for D2 forest
#' @param treeNo  tree No, need to be >= 50, if =10 has this error: Error in lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,  : 0 (non-NA) cases
#' @param iterationNo  iteration No
#' @param tune_unit tune unit  
#' 
#' @return the tuned leaf size 
#'
#' @export
 
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
    
    iCF_D2_exp = iCF(leafsize_initial$D2,  treeNo, tunepara = "none", iterationNo, Ntrain, "D2",  split_val_round_posi)
    iCF_D2_BT_exp <<- GET_TREE_L(iCF_D2_exp, 1) 
  }

  
  
  iCF_D3_exp = iCF(leafsize_initial$D3,  treeNo, tunepara = "none", iterationNo, Ntrain, "D3",  split_val_round_posi)
  iCF_D3_BT_exp <<- GET_TREE_L(iCF_D3_exp, 1) 
  
  iCF_D4_exp = iCF(leafsize_initial$D4,  treeNo, tunepara = "none", iterationNo, Ntrain, "D4",  split_val_round_posi)
  iCF_D4_BT_exp <<- GET_TREE_L(iCF_D4_exp, 1) 

  iCF_D5_exp = iCF(leafsize_initial$D5,  treeNo, tunepara = "none", iterationNo, Ntrain, "D5",  split_val_round_posi)
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





#'CF_RAW_key
#'
#' This function is the wrapper function to run raw causal forest
#' @param Train_cf training data
#' @param min.split.var minimum number of variables to split population
#' @param variable_type high-dimensional (HD) or non-HD
#' @param hdpct if high-dimensional (HD) variable, the top percentile used to split population
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
  #28=7*4, i.e., 7 is the minimum node# for a depth 4 Symmetrical tree, and when 7 variables are the q4 part of distribution, it requires 28 varaibles  
  #60=15*4, i.e., 15 is the minimum node# for a depth 5 Symmetrical tree, and when 7 variables are the q4 part of distribution, it requires 28 varaibles 
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


#'GET_TREE_DF
#'
#' This function is the wrapper function to run raw causal forest
#' @param No_nodes number of nodes from a best tree from a causal forest
#' @param treetype list format of a tree from causal forest
#' @param split_val_round_posi split value rounding position
#' 
#' @return the key results from raw causal forest 
#'
#' @export

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






#' GET_ORI_TREE_L
#' 
#' This function connects the majority vote (or any) deicison with the list of original tree and ouptut its list format
#' @param Maj_Tree_Decision a majority voted tree decision df that need to be presented in a tree plot 
#' @param ori_TreeL  #a list of original tree   
#' 
#' @return The original tree (lsit) that match the voted tree decision (df) 
#'
#' @export

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

#' HTE_P_extract
#' 
#' This function connects the majority vote (or any) deicison with the list of original tree and ouptut its list format
#' @param iCF_D output of iCF at a certain depth 
#' @param tree_depth  "D2", "D3", "D4", or "D5"   
#' 
#' @return P_value for heterogeneity of each causal forest 
#'
#' @export
 
HTE_P_extract <- function(iCF_D, tree_depth){
  HTE_P_cf_list <- lapply(iCF_D, function(df) mean(df$HTE_P_cf))
  HTE_P_cf   <- as.data.frame(rlist::list.rbind(HTE_P_cf_list)) #Bind all list elements by row
  colnames(HTE_P_cf) <- c("HTE_P_cf")
  HTE_P_cf$tree_depth <- tree_depth
  HTE_P_cf <- tibble::rowid_to_column(HTE_P_cf, "k")
  return(HTE_P_cf)
}
