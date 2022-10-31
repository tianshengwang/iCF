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
GRID_LEAFSIZE <- function (depth, initial, iCF_D_BT_exp_SG_D_mean, tune_unit){
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
    grid_MinNodeSize_D <<- seq(initial, initial*4, by=tune_unit*2)
  } else if ( iCF_D_BT_exp_SG_D_mean > depth #need a shallower tree
              #=================
              #Scenario 3: suggested value is too large
              #minimum leaf size for D5 (N/200) is too small to grow a very shallows, D5 tree,need to have a larger (smaller than initial value 200) minimum leaf size for D5, thus will count from 200 down to 100
              #=================
  ) {
    grid_MinNodeSize_D <<- seq(initial, initial/2, by=-tune_unit)
  }
  
  #########################################
  # STEP 2: try each sample size one by one
  #########################################
  for (m in 1: length(grid_MinNodeSize_D) ){
    iCF_D_experiment <- iCF(grid_MinNodeSize_D[m],  treeNo, iterationNo, Ntrain, paste0("D",depth), 1)
    iCF_D_BT_experiment <- GET_TREE_L(iCF_D_experiment, 1) 
    
    iCF_D_BT_experiment_SG        <- lapply(iCF_D_BT_experiment, function(df)   TREE2SUBGROUP(df)$majority)
    iCF_D_BT_experiment_SG_length <- lapply(iCF_D_BT_experiment_SG, function(df)   as.data.frame(df) %>% dplyr::mutate(length=stringr::str_count(subgroup, "&")  + 2 ) )
    iCF_D_BT_experiment_SG_Depth  <- lapply(iCF_D_BT_experiment_SG_length, function(df)   as.data.frame(df) %>% dplyr::mutate(Depth=max(length)   ) )
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
LEAFSIZE_tune <- function(Ntrain, initial_D2, treeNo, iterationNo){
  
  #suggested N to start with, then check, if doesn't work then tune by increasing or  minimum leafsize
  leafsize_initial = list(D5=initial_D2*5,
                          D4=initial_D2*4, 
                          D3=initial_D2*2,
                          D2=initial_D2)
  
  #=================
  #CHECK: 
  #Scenario 1: suggested value works! if node # of D2 tree==3, then return leafsize, otherwise rerun tuning!
  #=================
  #########################
  # Checking
  #########################
  iCF_D2_exp = iCF(leafsize_initial$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi)
  iCF_D2_BT_exp <<- GET_TREE_L(iCF_D2_exp, 1) 
  
  iCF_D3_exp = iCF(leafsize_initial$D3,  treeNo, iterationNo, Ntrain, "D3",  split_val_round_posi)
  iCF_D3_BT_exp <<- GET_TREE_L(iCF_D3_exp, 1) 
  
  iCF_D4_exp = iCF(leafsize_initial$D4,  treeNo, iterationNo, Ntrain, "D3",  split_val_round_posi)
  iCF_D4_BT_exp <<- GET_TREE_L(iCF_D4_exp, 1) 

  iCF_D5_exp = iCF(leafsize_initial$D5,  treeNo, iterationNo, Ntrain, "D3",  split_val_round_posi)
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
  
  
  D5_tune <- GRID_LEAFSIZE(5, initial_D2*5, iCF_D5_BT_exp_SG_D_mean, 10)
  D4_tune <- GRID_LEAFSIZE(4, initial_D2*4, iCF_D4_BT_exp_SG_D_mean, 10)
  D3_tune <- GRID_LEAFSIZE(3, initial_D2*2, iCF_D3_BT_exp_SG_D_mean, 10)
  D2_tune <- GRID_LEAFSIZE(2, initial_D2,   iCF_D2_BT_exp_SG_D_mean, 10)
  
    #########################
    # STEP 3. double check
    #########################
    leafsize_tune <- list(D5=D5_tune, #takes no computation time
                          D4=D4_tune, #takes no computation time
                          D3=D3_tune, 
                          D2=D2_tune
    )
   return(leafsize_tune)
}
