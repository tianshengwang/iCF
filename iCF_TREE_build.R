#' This function connects the majority vote (or any) deicison with the list of original tree in list format
#' @param Maj_Tree_Decision a majority voted tree decision df that need to be presented in a tree plot 
#' @param ori_TreeL  #a list of original tree list   
#' 
#' @return The original tree (lsit) that match the voted tree decision (df) 
#'
#' @export
#' 

GET_TREE_L <- function (iCF_D, listindex) {
  iCF_D <- lapply(iCF_D,"[", listindex) 		#subsetting nested list to get besttree, extract the "listindex"th list from mother list
  if (listindex == 1 | listindex== 3 ){     # listindex=1 to extract besttreeelist, listindex=3 to extract splitfreqlist
    iCF_D_2 <- lapply(iCF_D, function(L) data.frame(L[[1]]) ) 
  } else {
    iCF_D_2 <- lapply(iCF_D, function(L) L[[1]] ) 
  }
  return(iCF_D_2)
}



#' This function tune leaf size for iCF 
#' @param Ntrain sample size
#' @param treeNo  tree No
#' @param iterationNo  iteration No  
#' 
#' @return the tuned leaf size 
#'
#' @export
#' 
LEAFSIZE_tune <- function(Ntrain, treeNo, iterationNo){
  
    #suggested N to start with, then check, if doesn't work then tune by increasing or  minimum leafsize
    leafsize = list(D4=100, 
                    D3=50,
                    D2=25)
    
    MeanNodeNum_D2_L <- list()
    MeanNodeNum_D3_L <- list()
    #=================
    #CHECK: 
    #Scenario 1: suggested value works! if node # of D2 tree==3, then return leafsize, otherwise rerun tuning!
    #=================
    #########################
    # Checking
    #########################
    iCF_D2_exp = iCF(leafsize$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi)
    iCF_D2_BT_exp <<- GET_TREE_L(iCF_D2_exp, 1) 
    if ( round(mean( do.call("rbind" , lapply(iCF_D2_BT_exp, function(df) nrow(df))) ), 0) ==3){
      return(leafsize)
    
    } else {
      
  #########################################
  # STEP 1: get a list of suggested value
  #########################################
      if ( round(mean( do.call("rbind" , lapply(iCF_D2_BT_exp, function(df) nrow(df))) ), 0) < 3 
    #=================
    #Scenario 2: suggested value is too small
    #minimum leaf size for D2 (N/25) is too large to split,need to have a smaller (larger than initial value 25) minimum leaf size for D2, thus will count from 25 up to 100
    #=================
                ) {
    grid_MinNodeSize_D2 <<- seq(25, 100, by=5)
    } else if ( round(mean( do.call("rbind" , lapply(iCF_D2_BT_exp, function(df) nrow(df))) ), 0) >3 
    #=================
    #Scenario 3: suggested value is too large
    #minimum leaf size for D2 (N/25) is too small to grow a very shallows, D2 tree,need to have a larger (smaller than initial value 25) minimum leaf size for D2, thus will count from 25 down to 1
    #=================
    ) {
    grid_MinNodeSize_D2 <<- seq(25, 1, by=-5)
    }
    #########################################
    # STEP 2: try each sample size one by one
    #########################################
      for (m in 1: length(grid_MinNodeSize_D2) ){
        iCF_D2_experiment <- iCF(grid_MinNodeSize_D2[m],  200, 50, Ntrain, "D2", 1)
        iCF_D2_BT_experiment <- GET_TREE_L(iCF_D2_experiment, 1) 
        MeanNodeNum_D2_L[[m]] <-round(mean( do.call("rbind" , lapply(iCF_D2_BT_experiment, function(df) nrow(df))) ), 0)
        
        iCF_D3_experiment <- iCF(grid_MinNodeSize_D2[m]*2,  200, 50, Ntrain, "D3", 1)
        iCF_D3_BT_experiment <- GET_TREE_L(iCF_D3_experiment, 1) 
        MeanNodeNum_D3_L[[m]] <-round(mean( do.call("rbind" , lapply(iCF_D3_BT_experiment, function(df) nrow(df))) ), 0)
        
        if (MeanNodeNum_D2_L[[m]] %in% c(3) & 
            MeanNodeNum_D3_L[[m]] %in% c(7)  )     break #if mean node # of best trees from iCF equal to 3, i.e. D2 tree structure
        
      }
    #########################
    # STEP 3. double check
    #########################
      leafsize <- list(D4=grid_MinNodeSize_D2[length(MeanNodeNum_D2_L)]*4, 
                      D3=grid_MinNodeSize_D2[length(MeanNodeNum_D2_L)]*2, 
                      D2=grid_MinNodeSize_D2[length(MeanNodeNum_D2_L)]
                      )
      #testing for tuned N 
      iCF_D2_exp = iCF(leafsize$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi)
      iCF_D2_BT_exp <<- GET_TREE_L(iCF_D2_exp, 1) 
      if (round(mean( do.call("rbind" , lapply(iCF_D2_BT_exp, function(df) nrow(df))) ), 0)==3){
        return(leafsize)
      }  else {
        stop("Automatic tuning minimum leaf size is not successful, need manual adjustment")
      }
   #---------- double check over -----------------------

    }
}






CF_RAW_key <- function(Train_cf){
  s_rawCF =  Sys.time()
  X <- Train_cf[,vars_forest]
  Y <- as.vector( as.numeric( Train_cf[,1] ) )
  W <- as.vector( as.numeric( Train_cf[,2] ) )  
  #------------------------------------------------------------------------
  # directly use the following for EMPIRIAL STUDIES, i.e. Real-World Data
  #------------------------------------------------------------------------
  #------------------------------------------------------------------------
  W.forest <- regression_forest(X, W)
  W.hat    <- predict(W.forest)$predictions
  
  #forest-based method to predict Y
  Y.forest <- regression_forest(X, Y)
  Y.hat    <- predict(Y.forest)$predictions
  
  cf.raw = causal_forest(X, Y,  W,  Y.hat = Y.hat, W.hat = W.hat)
  
  best_tree_info<-find_best_tree(cf.raw, "causal")
  best_tree_info$best_tree
  # Plot trees
  par(mar=c(1,1,1,1))
  tree.plot = plot(grf::get_tree(cf.raw, best_tree_info$best_tree))
  tree.plot
  
  
  #--------------------------------------------------------------
  #get P-value of omnibus test for HTE presence 
  #---------------------------------------------------------------
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
  selected_cf.idx_median = which(varimp_cf > quantile(varimp_cf, 0.5) )
  selected_cf.idx_q3 = which(varimp_cf > quantile(varimp_cf, 0.75) )
  
  #28=7*4, i.e., 7 is the minimum node# for a depth 4symmetrical tree, and when 7 variables are the q4 part of distribution, it requires 28 varaibles  
  if (length(selected_cf.idx_mean)==1 & 
      dim(X)[2]>=28# number of covariates >= 28
  ) {
    selected_cf.idx = selected_cf.idx_q3
  } else if( length(selected_cf.idx_mean)==1 & 
             dim(X)[2] < 28 # number of covariates > 28
  ) {
    selected_cf.idx = selected_cf.idx_median
  } else if( length(selected_cf.idx_mean) > 1 ) {
    selected_cf.idx = selected_cf.idx_mean
  }
  #---------------------------------------------------------------
  head( X[, selected_cf.idx], 1 ) #be careful!!! not head(Train[, selected_cf.idx])
  
  colnames(X[, selected_cf.idx ] )
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
