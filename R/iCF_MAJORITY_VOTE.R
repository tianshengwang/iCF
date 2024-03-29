#############################################################################################
# MAJORITY VOTE
# Author: Tiansheng Wang  
# Last update date:7/29/2023
# Version: 1.0         
#############################################################################################

#################################################----------------------------------------------
#################################################----------------------------------------------
# IV. majority vote
#################################################----------------------------------------------
#################################################----------------------------------------------
#
#' MAJORITY_COUNT
#' 
#' This function count the occurence of majority 
#' @param list_ref list of unique trees, tree relaxed, or subgroups
#' @param list_all prepared list of trees, tree relaxed, or subgroups
#' 
#' @return The occurency of majority
#'
#' @export
#'
MAJORITY_COUNT <- function(list_ref, list_all){
  n_occurence <- NULL
  for (elt in list_ref){
    n_occ <- 0 # A counter of occurence for this e
    for (df in list_all){
      if (identical(elt, df)){ # If df is exactly elt, we count it
        n_occ <- n_occ + 1
      }
    }
    n_occurence <- c(n_occurence, n_occ)
  }
  return(n_occurence)
}

#' STANDARD_CHECK
#' 
#' This function connects if identical or not label to a list and subset the list to exclude those not identical to reference list
#' if_ID means if Identical to MAJORITY VOTED or TRUTH
#' @param standardLIST a label list of being identical or not 
#' @param inputLIST the list need to be labeled and subsetted 
#' 
#' @return The subsetted list that all are identical to reference list.
#'
#' @export

STANDARD_CHECK <- function(standardLIST, inputLIST){
  if_ID_2_vote  <- rlist::list.zip(identical=standardLIST, df=inputLIST)  
  ID_2_vote     <- rlist::list.exclude(if_ID_2_vote, identical=='FALSE')    #exclude dataframes not identical with majority vote
  pure_identical<-lapply(ID_2_vote , function(df) rlist::list.extract(df, 2)) #remove "TRUE" logic item from lists in the list
  return(pure_identical)
}


#' MAJORITY_PICK_LEAF
#' 
#' This function pick those match matjority vote and with LEAF or NODE
#' @param inputlist a list of df match the voted
#' @param leaf_para "TRUE", means it is a leaf
#' 
#' @return The subsetted list that match voted
#'
#' @export
#'
#' @example:  all_identical_leaf <- MAJORITY_PICK_LEAF(if_ID_Y_W, "TRUE")

MAJORITY_PICK_LEAF <- function (inputlist, leaf_para) {
  pure_identical_LorN <- lapply(inputlist, function(df) subset(df, is_leaf== leaf_para)) #extract LEAF node
  all_identical_LorN  <- rlist::list.stack(pure_identical_LorN) #combine all dataframes in the list into a dataframe, make a data frame of voted majority tree's NODE
  return(all_identical_LorN)
}

#' SG_PROBABILITY
#' 
#' This function check the probability of the majority voted subgroup from iCF 
#' @param list1 prepared list of trees, tree relaxed, or subgroups
#' 
#' @return The subsetted list that all are identical to reference list.
#'
#' @export

SG_PROBABILITY <- function( list1){
  #get a list of unique Tree (dataframe), Identify the unique dataframe
  list1_u <- unique(list1) #this function get unique SG
  N_occur_majority <- MAJORITY_COUNT(list1_u, list1) #This function count the occurence of majority 
  N_occur_prob <- N_occur_majority/sum(N_occur_majority)
  return(N_occur_prob)
}


#' MAJORITY_VOTE
#' 
#' This function performs majority vote among best trees from iterative CF 
#' @param list0 original list of trees
#' @param list1 prepared list of trees, tree relaxed, or subgroups
#' @param list2 only useful for prepared list if by subgroup, then remove node column to count for majority.tree
#' @param list3 original list of split frequency!
#' @param truth truth 
#' @param truth_N1 true top node for tree
#' @param truth_N123 true tope 3 nodes for tree
#' @param split_val_round_posi rounding parameter for decimal place
#' 
#' @return The subsetted list that all are identical to reference list.
#'
#' @export

MAJORITY_VOTE <- function(list0, #original list of trees
                          list1, #prepared list of trees, tree relaxed, or subgroups
                          list2, #only useful for prepared list if by subgroup, then remove node column to count for majority.tree
                          list3, #original list of split frequency!
                          truth, 
                          truth_N1, #for tree
                          truth_N123, #for tree
                          split_val_round_posi
){
  #----------------------------------------------------------------------------
  #STEP1: compared to the prepared list, i.e., list1, to get majority
  #----------------------------------------------------------------------------		  
  
  #first things first: get ncol to automatically distinguish among tree, tree_r(ncol=5), or subgroup (ncol=2) input
  ncol_list1 <- lapply(list1, function(df) ncol(df))
  #round the mean value of columns
  mean_ncol_list1 <- round( mean(unlist(ncol_list1)) )
  #get a list of unique Tree (dataframe), Identify the unique dataframe
  list1_u <- unique(list1) 
  #count the occurence of majority 
  
  
  #get frequency of tree structure occurence
  N_occur_majority <- MAJORITY_COUNT(list1_u, list1)
  #get majority tree structure
  majority.tree <- list1_u[[which.max(N_occur_majority)]]
  #get majority tree structure, same as above
  majority.tree.1st <- list1_u[[  match(Rfast::nth(N_occur_majority, 1, descending = T), N_occur_majority  )]]
  #get 2nd popular tree structure 
  majority.tree.2nd <- list1_u[[   match(Rfast::nth(N_occur_majority, 2, descending = T), N_occur_majority  )]]
  #get 3rd popular tree structure
  majority.tree.3rd <- list1_u[[   match(Rfast::nth(N_occur_majority, 3, descending = T), N_occur_majority  )]]  
  #make truth a list (tree_true or tree_true_r)
  truth_list      <- list(truth)
  #get number of occurrence of truth in list1 (the prepared list of trees, tree relaxed, or subgroups)
  #N_occur_truth could be 0, e.g. truth = D4 tree, majority is D3 tree.
  N_occur_truth   <-  MAJORITY_COUNT(truth_list, list1)
  #make the top (root) node a list
  truth_N1_list   <- list(truth_N1)   
  #make the top 3 nodes a list, if truth=D2 tree, then "NA" as there will be only ONE node.
  truth_N123_list <- list(truth_N123)
  #if tree_r, then need to additionally drop split value
  if (mean_ncol_list1==5) { #i.e. if tree_r which has 5 columns: 
    list0_N1   <-lapply(list0, function(df) subset(df, node=="node-01", select = c(1:5) ) )
    list0_N123 <-lapply(list0, function(df) subset(df, node=="node-01"|node=="node-02"|node=="node-03", select = c(1:5) ) )
    
  } else {
    list0_N1   <-lapply(list0, function(df) subset(df, node=="node-01", select = c(1:6)))
    list0_N123 <-lapply(list0, function(df) subset(df, node=="node-01"|node=="node-02"|node=="node-03", select = c(1:6)))  

  }
  #get occurence of node-1 in the truth
  N_occur_truth_N1   <-  MAJORITY_COUNT(truth_N1_list,   list0_N1)
  #get occurence of node-123 in the truth, if truth_N123_list="NA", then this number will be ZERO.
  N_occur_truth_N123 <-  MAJORITY_COUNT(truth_N123_list, list0_N123)
  #-----------------------------------------------------------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------------------------------------------------------
  #check the voted majority tree is equal to truth or not
  majority_EQ_true      <- identical(majority.tree, truth)                #IDENTICAL also works, but ALL.EQUAL provides reasons when  it is FALSE
  #the following explain why if false:
  majority_EQ_true_full <- all.equal(majority.tree, truth)                
  #check if each single best tree(or subgroup) == the TRUE tree(subgroup)
  if_true     <-lapply(list1, function(df) identical(df, truth))          
  #check if each single best tree(or subgroup) == the majoriy tree(subgroup)
  if_majority <-lapply(list1, function(df) identical(df, majority.tree))  
  
  #------------------------------------------------------------------------------------------
  #STEP2: merge with orginal list of SPLITFREQUENCY, to get split frequentcy
  #------------------------------------------------------------------------------------------		  
  #save a list of split frequency for those identical to majority vote
  SF_ID_list <- STANDARD_CHECK(if_majority, #if_majority provides info on if df is identical to majority voted tree/tree_r/subgroup
                               list3)
  #------------------------------------------------------------------------------------------
  #STEP3: merge with orginal list of BESTTREE, i.e., list0, to get average hat{Y}, hat{W}
  #------------------------------------------------------------------------------------------		  
  #---------------------------------------------------------------------------------------------------------------------------------
  # step 3.1: TELL if running subgroup or not, if so merge with subgroup list to reorder node!
  if (mean_ncol_list1==2) {#i.e., running subgroup
    list0_subgroup   <- lapply(seq_along(list0), function(df) dplyr::left_join(list0[[df]], #list0: original list of trees
                                                                               list2[[df]], #list2: only useful for prepared list if by subgroup,
                                                                               by=c("node"="LEAF")))	
    #To avoid specific column name (k for iCF, b for oneCF):
    #list0_subgroupID <- lapply(list0_subgroup,   function(df) subset(df, select=c(subgroupID, is_leaf, left_child, right_child, split_variable, split_value, samples, avg_Y, avg_W, HTE_P_cf, k)) )
    list0_subgroupID <- lapply(list0_subgroup,   function(df) subset(df, select=c(12, 2:11)) )
    
    #need to assess average of hat{Y}, hat{W} and average of split_value spearately for subgroup
    #the former can be assessed in subgroup, the latter (i.e., split_value) can only be assessed in tree list
    if_ID_Y_W       <- STANDARD_CHECK(if_majority, list0_subgroupID)  #give name to each list, e.g. name if_majority as "identical"
    if_ID_split_val <- STANDARD_CHECK(if_majority, list0) 
  } else {
    if_ID_Y_W       <- STANDARD_CHECK(if_majority, list0)
    if_ID_split_val <- if_ID_Y_W
  }
  
  
  #-----------------------------------------------------------------------------------------------------------------------------		  
  #make a data frame of all those df that matches voted majority tree's LEAF
  all_identical_leaf <- MAJORITY_PICK_LEAF(if_ID_Y_W, "TRUE")
  #convert to numeric
  all_identical_leaf[c("samples", "avg_Y", "avg_W")] <- sapply(all_identical_leaf[c("samples", "avg_Y", "avg_W")],as.numeric) 
  
  #drop character variable and obtain average value of sample size, Y, and W for each leaf
  #group_by doesnt' work possibly due the the way I generated data from list!!! wrong code: all_identical_leaf_small %>% group_by(node) %>% summarise_at(vars("samples"), mean)
  if (mean_ncol_list1==2) { #i.e. if subgroup list
    avg_leaf <- aggregate(cbind(samples, avg_Y, avg_W) ~ subgroupID, data = all_identical_leaf, mean, na.rm = TRUE) #rounding doesn't work, give up here, do it separately
  } else {
    avg_leaf <- aggregate(cbind(samples, avg_Y, avg_W) ~ node,       data = all_identical_leaf, mean, na.rm = TRUE) #rounding doesn't work, give up here, do it separately
  }
  #rounding:
  avg_leaf$samples<-round(avg_leaf$samples,0)
  avg_leaf$avg_Y  <-round(avg_leaf$avg_Y,2)
  avg_leaf$avg_W  <-round(avg_leaf$avg_W,2)
  #-----------------------------------------------------------------------------------------------------------------------------
  #STEP4: merge with orginal list, i.e., list0, to get average split value. Have to rely on tree structure, even for subgroup
  #-----------------------------------------------------------------------------------------------------------------------------			  
  #make a data frame of all df of "voted majority" tree's NODE!
  all_identical_node <- MAJORITY_PICK_LEAF(if_ID_split_val, "FALSE")
  #convert to numeric
  all_identical_node[c("split_value")] <- sapply(all_identical_node[c("split_value")],as.numeric) 
  #----------------------------------------------------
  #plot split_value by node
  split_value_plot <-  ggplot2::ggplot(all_identical_node, #[which(all_identical_node$node !='node-1'),]
                                       ggplot2::aes(x = split_value, y = node)) + 
    ggridges::geom_density_ridges(ggplot2::aes(fill = node)) 
  # ggsave(paste("split_value_plot", minnodesize, treeNo, iterationNo, Ntrain, tree_depth, intTRUE, name, ".png", sep="_"))
  #----------------------------------------------------
  #get HTE P values from those trees have the same structure as majority vote
  # "k" for iCF, "b" for oneCF
  if         ( dplyr::nth(colnames(all_identical_node), 10) == "k") { 
    HTE_P_cf_M        <- aggregate(cbind(HTE_P_cf) ~ k,    data = all_identical_node, mean,     na.rm = TRUE)
    HTE_P_cf_M_median <- median(HTE_P_cf_M$HTE_P_cf) 
  } else if  ( dplyr::nth(colnames(all_identical_node), 10) == "b") {
    HTE_P_cf_M        <- aggregate(cbind(HTE_P_cf) ~ b,    data = all_identical_node, mean,     na.rm = TRUE)
    HTE_P_cf_M_median <- all_identical_node[1,11]
  }
  
  #get mean, SD, distribution: 
  split_value_mean     <- aggregate(cbind(split_value) ~ node, data = all_identical_node, mean,     na.rm = TRUE) %>% 
    dplyr::mutate(split_value_mean = round(split_value,  split_val_round_posi)) %>%
    dplyr::rename(split_value_mean_raw = split_value)
  
  split_value_sd       <- aggregate(cbind(split_value) ~ node, data = all_identical_node, sd,       na.rm = TRUE) %>%
    dplyr::mutate(split_value_sd     = round(split_value,    split_val_round_posi)) %>%
    dplyr::rename(split_value_sd_raw = split_value)
  
  split_value_quantile <- aggregate(cbind(split_value) ~ node, data = all_identical_node, quantile, na.rm = TRUE) %>%
    dplyr::mutate(across(where(is.numeric), round, split_val_round_posi))
  #won't calculate MODE as there might be multiple MODE number 

  #----------------------------------------------------------------------------------------
  #obtain a synthetic tree by combing tree structure (split_variable only) with mean split_value
  majority.syn0 <-  dplyr::left_join(majority.tree %>% select(c(1:5)), #make sure drop split_value from original df for tree (i.e., when not for tree_r)
                                     dplyr::select(split_value_mean,   #left join with df for  mean split 
                                                   c(node, split_value_mean)),
                                     by = "node" )  %>%
                    dplyr::rename(split_value=split_value_mean) %>%    #rename mean split value
                    dplyr::left_join(dplyr::select(avg_leaf,           #left joint df for sample size, avg_Y, avg_W
                                                   c(node, samples, avg_Y, avg_W)),
                                     by = "node" ) %>%
                    dplyr::mutate(k=1, 
                                  HTE_P_cf = mean(HTE_P_cf_M$HTE_P_cf))
  
  majority.syn1 <- majority.syn0 %>% dplyr::filter(split_variable != "NA") 
  
  #In case the synthetic mean split value is not integer for ordinal variable (bianry variable is perfectly fine as it always split at "0")
  # test code: assign a non integer to ordinal X2's split value, e.g., by: majority.syn1 [,"split_value"][which(majority.syn1$node == "node-01" )] <- 1.7; 
  # then round applied to each row of the sub dataset
  rounded_split_value <- apply(majority.syn1, 1, function(x) if ( length(unique (  eval(parse(text=paste("X$", x[5], sep = ""))) )) <8) {
    round(  as.numeric(x[6]), digits = 0) 
  } else {
    as.numeric( x[6] )
  }
  )
  #replace the unrounded split value with the round one
  majority.syn1$split_value =  rounded_split_value
  #pick key columns
  majority.syn2 <- majority.syn1 %>% dplyr::select(node, split_variable, split_value)
  #----------------------------------------------------------------------------------------
  
  
  #Replace split_value from majority.syn2 (synthetic tree's average split_value) by node 

  majority.syn0 <-  dplyr::left_join(majority.syn0, majority.syn2, by = c("node", "split_variable")) %>% 
                    dplyr::mutate(Values = ifelse(!is.na(split_value.x), split_value.y, split_value.x)) %>% 
                    dplyr::select(-c(split_value.x, split_value.y)) %>%
                    dplyr::rename(split_value=Values) %>%
                    dplyr::relocate(split_value, .after=split_variable)
  
  # devide dataframe into several dataframe by row and cobmine in a LIST
  #_______________________________________________
  stability     <- max(N_occur_majority)/sum(N_occur_majority)
  stability.2nd <- Rfast::nth(N_occur_majority, 2, descending = T)/sum(N_occur_majority)
  stability.3rd <- Rfast::nth(N_occur_majority, 3, descending = T)/sum(N_occur_majority)
  accuracy      <- N_occur_truth     /sum(N_occur_majority)
  accuracy_N1   <- N_occur_truth_N1  /sum(N_occur_majority)
  accuracy_N123 <- N_occur_truth_N123/sum(N_occur_majority)
  max_occurence <-max(N_occur_majority)
  
  majority.list <- list(n_occurence          =N_occur_majority, 
                        max_occurence        =max_occurence, 
                        stability            =stability,  #=max/total
                        stability.2nd        =stability.2nd,  #=2nd max/total
                        stability.3rd        =stability.3rd,  #=3rd max/total
                        accuracy             =accuracy,  #=accurate/total
                        accuracy_N1          =accuracy_N1,  #=accurate/total
                        accuracy_N123        =accuracy_N123,  #=accurate/total
                        majority_EQ_true     =majority_EQ_true, 
                        majority_EQ_true_full=majority_EQ_true_full, 
                        majority             =majority.tree, 
                        majority.2nd         =majority.tree.2nd, 
                        majority.3rd         =majority.tree.3rd, 
                        majority.syn         =as.data.frame(majority.syn0),
                        truth                =truth,
                        avg_leaf             =avg_leaf, 
                        split_value_mean     =split_value_mean, 
                        split_value_sd       =split_value_sd,  
                        split_value_quantile =split_value_quantile,
                        split_value_plot     =split_value_plot,
                        HTE_P_cf_M_median    =HTE_P_cf_M_median,
                        HTE_P_cf_M           =HTE_P_cf_M,
                        SF_ID_LIST           =SF_ID_list)
  return(majority.list)
}



#' CV_SG_MAJORITY
#' 
#' This function finds the majority voted subgroup decision across CV
#' @param vote_D_subgroup.L list of vote_D_subgroup across CV
#' @param stability_D_T_r stability of tree structures ignoring spliting value at each depth
#' @param K cross-validation fold 
#' 
#' @return The majority voted subgroup decision across CV
#'
#' @export

CV_SG_MAJORITY <- function(vote_D_subgroup.L, stability_D_T_r, K){
  
  #This function count the occurence of majority 
  N_occur_majority_CV <- MAJORITY_COUNT(unique(vote_D_subgroup.L),vote_D_subgroup.L)
  
  #if none of the SG across CV is the same, pick the one with the largest vote stability score
  if(length(N_occur_majority_CV )==K ){
    selected_CV_idx=which.max(stability_D_T_r) #prevent randomly selecting the 1st one in the CV list
    majority_CV=vote_D_subgroup.L[[selected_CV_idx]]
  } else {
    #find the selected
    selected_CV_idx_pre = which.max(N_occur_majority_CV)
    majority_CV <- unique(vote_D_subgroup.L)[[ selected_CV_idx_pre ]]
    
    selected_CV_idx = which(lapply(vote_D_subgroup.L, function(x) identical(majority_CV, x) )==TRUE) #will be mutiple, i.e.>=2
    
  }
  
    stability_CV     <- max(N_occur_majority_CV)/sum(N_occur_majority_CV)

    return(list(majority_CV=majority_CV,
              selected_CV_idx=selected_CV_idx,
              stability_CV=stability_CV))
}


#' CVBIAS_D_MAJORITY
#' 
#' This function returns the mean MSE value of those majority voted subgroup decisions across CV
#' @param vote_D_subgroup.L list of voted subgroup decisions, of those not majorityed, their MSE will be ignored
#' @param stability_D_T_r stability of tree structures ignoring spliting value at each depth
#' @param K cross-validation fold 
#' @param metric list of metric including RMSE, Rsquared, MAE for continuous variables
#' 
#' @return The the mean MSE value across CV
#'
#' @export
#'
#'
CVBIAS_D_MAJORITY <- function(vote_D_subgroup.L, stability_D_T_r, K, metric){
  #if probability of majority voted structure is 1
  #SG_PROBABILITY: check the probability of the majority voted subgroup from iCF 
  if (length(SG_PROBABILITY (vote_D_subgroup.L) ) ==1 ){
    mse_all_folds <- as.data.frame(do.call("rbind", metric))
    majority.SG = vote_D_subgroup.L[[1]]
  } else {
    
    if_majority <-lapply(vote_D_subgroup.L, 
                         function(df) identical(df, CV_SG_MAJORITY(vote_D_subgroup.L,
                                                                   stability_D_T_r, 
                                                                   K)$majority_CV #finds the majority voted subgroup decision across CV 
                                                                   ))  #if none of SG appear > once, then the 1st one will be picked as the "majority"
    #connects if identical or not label to a list and subset the list to exclude those not identical to reference list
    #if none of SG appear > once, then the 1st one will be picked 
    mse_per_fold_puremajor <- STANDARD_CHECK(if_majority,metric)
    mse_all_folds <- as.data.frame(do.call("rbind",mse_per_fold_puremajor))
    
  }
  #----------------------------------------------------------------------------------
  mse_all_folds <- cbind(mse_all_folds, "MSE"=(mse_all_folds[,"RMSE"])^2) 
  mse_all_folds <- mse_all_folds %>%
    dplyr::mutate("CV_MSE"=mean(mse_all_folds[,"MSE"]),
                  "CV_MSE_SE"=sd(mse_all_folds[,"MSE"]))
  
  return(mse_all_folds$CV_MSE[1])
}

#' CVBIAS_MAIN
#' 
#' This function returns the mean MSE value across CV
#' @param metric list of metric including RMSE, Rsquared, MAE for continuous variables
#' 
#' @return The the mean MSE value across CV
#'
#' @export
#'

CVBIAS_MAIN <- function(metric){
  mse_all_folds <- as.data.frame(do.call("rbind", metric))
  mse_all_folds <- cbind(mse_all_folds, "MSE"=(mse_all_folds[,"RMSE"])^2) 
  mse_all_folds <- mse_all_folds %>%
    dplyr::mutate("CV_MSE"=mean(mse_all_folds[,"MSE"]),
                  "CV_MSE_SE"=sd(mse_all_folds[,"MSE"]))
  
  return(mse_all_folds$CV_MSE[1]) 
}
