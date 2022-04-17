#############################################################################################
# Parent Node
# Author: Tiansheng Wang  
# Last update date:12/20/2020
# Version: 0.1         
#############################################################################################
######################################################################################################----------------------------------------------
######################################################################################################----------------------------------------------
########                      III. prepare subpopulation-based majority vote. PART I.
######################################################################################################----------------------------------------------
######################################################################################################-----------------------------------------------
#' @param Child_node children nodes of upper level nodes (nodes at a shallower depth)
#' @param tree_original the list of best trees from iCF (iCF_D4_BT, iCF_D3_BT, iCF_D2_BT)
#' @param tree the ouput of previous parent level, if parent_level=1 then tree=tree_original 
#' @param parent_level the level of parent node. The key idea of this function is backward stagewise method to obtain parent nodes of a upper level (shallower level) until reach node-01 
#' #' #---------------------------------------------------------------------------------------------------------------
# III. "PARENT" function must run independently, recursivley to get the parent node, so that we can get subgroup decision
# this is the most imporant and complicated part of the program: subgroup identification based on trees
# the basic idea is using child node# to identify parent node# recursively, i.e.up to 3 times as we only deal with 4-way interactions using depth=4 trees (the 4th level is for leaves)
#---------------------------------------------------------------------------------------------------------------
#lab start
#tree_original<-iCF_D2_BT
#tree         <-iCF_D2_BT
#parent_level=1
#-------------------------------
#tree_original<-iCF_D4_BT
#tree<-parent_final
#parent_level=2


#...
#lab over

#this PARENT_NODES is an intermediate function for PRE_MAJORITY_SUBGROUP
PARENT_NODES <- function(tree_original, tree, parent_level){
######################################################################################################################################################################  
######################################################################################################################################################################  
#SUBSECTION 1: [ORIGINAL] leaf info, depends on depth level thus included this part in the Fx,  from ORIGINAL tree, this part is constant for a specfic depth level!
######################################################################################################################################################################  
######################################################################################################################################################################  
#get all splitter and exclude leaves                                                                                                                                 #
splitter_L_original <- lapply(tree_original, function(df) df[which(df$is_leaf=='FALSE'), ] )                                                                         # 
#keep node info and drop other variables for splitters                                                                                                               #
splitter_node_L_ori <- lapply(splitter_L_original, function(df) df[c("node")] )                                                                                      #
#get all leaves and exclude splitters                                                                                                                                #
leaf_L_original <- lapply(tree_original, function(df) df[which(df$is_leaf=='TRUE'), ])                                                                               #                     
#keep node info and drop other variables for leaves                                                                                                                  #
leaf_node_L_ori <- lapply(leaf_L_original, function(df) df[c("node")] )                                                                                              #
#rename "node" as "LEAF"                                                                                                                                             #
LEAF_L_ori  <- lapply(leaf_node_L_ori, function(df) {names(df)[names(df) == 'node'] <- 'LEAF';df})                                                                   #
######################################################################################################################################################################  
######################################################################################################################################################################  
#SUBSECTION 2: [EXTRACT] the node # which we'll use to identify parent nodes	  
######################################################################################################################################################################  
######################################################################################################################################################################  
if (parent_level==1){
  #GET LEAFE NODE OF ORIGINAL INPUT----------------------------------------------------------------------------------------------------------------------------------#
  #keep node info and drop other variables read from leaf_list depends on parent_level                                                                               #
  leaf_node_list <- leaf_node_L_ori      #lapply(leaf_list,         function(df) df[c("node")] )                                                                     #
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------#
} else {
  #############################################################################################################################################################################################
  #Step 3.1 subsetting row 1 for "node-01"                                                                                                                                                    #
  #############################################################################################################################################################################################
  #extract the first row of parent_level, could be node=01 or empty (i.e. node-01 does not exist as the parent notde when parent level=1)                                                                                                                                                      #
  leaf_L_N1   <- lapply(tree, function(df)   df[which(df[,1] =="node-01" | df[,1] =="node-00001" ), ]   )                                                                                                                     #  
  #extract the non-Node-01 row as label to check if continue run or output                                                                                                                          #
  leaf_L_woN1 <- lapply(tree, function(df)   df[which(df[,1] !="node-01" & df[,1] !="node-00001" ), ]   )                                                                                                                   # 
  #extract actual "parent level" info
  leaf_L_N1_colname <- lapply(leaf_L_N1, function (df) as.numeric( stringr::str_sub(colnames(df[1]),-1,-1) ) )
  #rename, i.e. remove the suffix for "parent_level" version                                                                                                                                        #
  leaf_L_N1_rename_pre <- list.zip(leaf_L_N1_colname, leaf_L_N1)
  leaf_L_N1_rename <- lapply( leaf_L_N1_rename_pre, function(df) {names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('parent_left_child',  df$leaf_L_N1_colname, sep="_")] <-  'parent_left_child' ;                                        #   
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('parent_right_child', df$leaf_L_N1_colname, sep="_")] <-  'parent_right_child';                                        # 
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('parent_sign',        df$leaf_L_N1_colname, sep="_")] <-  'parent_sign'       ;                                        #             
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('parent_split_var',   df$leaf_L_N1_colname, sep="_")] <-  'split_variable'    ;                                        # 
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('parent_split_val',   df$leaf_L_N1_colname, sep="_")] <-  'split_value'       ;                                        # 
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('node',               df$leaf_L_N1_colname, sep="_")] <- 'node'               ;                                        #
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('Child_node',         df$leaf_L_N1_colname, sep="_")] <- 'Child_node'         ;                                        #
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('Freq',               df$leaf_L_N1_colname, sep="_")] <- 'Freq'               ;                                        #
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('left_child',         df$leaf_L_N1_colname, sep="_")] <- 'left_child'         ;                                        #
                                                                  names(df$leaf_L_N1)[names(df$leaf_L_N1) == paste('right_child',        df$leaf_L_N1_colname, sep="_")] <- 'right_child'        ; df$leaf_L_N1})                                   #  
  leaf_L_N1_ori <- lapply(tree_original,  function(df)   df[which(df[,1] =="node-01" | df[,1] =="node-00001" ), ]   )                                                                                                       # 
  leaf_L_N1_pre <- list.zip(leaf_L_N1_ori, leaf_L_N1_rename)                                                                                                                                          # 
 require(dplyr)
   leaf_L_N1_co <- lapply(leaf_L_N1_pre, 
                         function(df)  
                           dplyr::left_join(df$leaf_L_N1_ori, dplyr::select(df$leaf_L_N1_rename, c(node, parent_sign, parent_left_child, parent_right_child, Child_node)),by = "node" ) ) # 
  #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  #keep node (which is the parent node of previous output's Child_node) info and drop other variables read from leaf_list depends on parent_level
  leaf_node_list     <- lapply(tree,         function(df) df[1] )
  
  }
######################################################################################################################################################################  
######################################################################################################################################################################  
#SUBSECTION 3: [IDENTIFY] parent node # by the node number of of this current "parent level"	
######################################################################################################################################################################  
######################################################################################################################################################################  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#  
#step 3.0 prepare for key node# used for identificying parent node:
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#  
#extract frequency of node from parent_node
leaf_node_freq      <- lapply(leaf_node_list,     function(df) as.data.frame( table(df ) ) )                     
#EXCLUDE "node-01" as it cannot appear in left_child/right child column (making the same row # as )
leaf_node_freq_woN1 <- lapply(leaf_node_freq,    function(df) df [which(df$df!="node-01" & df$df!="node-00001"),] )    
#EXCLUDE "node-01" as it cannot appear in left_child/right child column (making the same row # as )
N_freq_woN1_nrow  <- lapply(leaf_node_freq_woN1,    function(df)  nrow(df)  )   
#Extract node number
leaf_node_No_L    <- lapply(leaf_node_list,     function(df) df$node_No <- sub('.*-', '', df$node))                                                                    # 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#  
#step 3.1 extract the node # of each leaf saving as a vector in a list                                                                                                 #  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------# 
#combine the node # of each leaf and splitter info                                                                                                                     #
List0 <- list.zip(leaf_node_No_L, splitter_L_original, N_freq_woN1_nrow, tree)                                                                                         #
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#  
#step 3.2 find the parent node by reading current "leaf_node_No_L" from left_child and right_child columns in original splitter info,respectively                      #
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------# 
List_left  <- lapply(List0, 
                     function(df) 
                     if (df$N_freq_woN1_nrow >0) {
                                                  df$splitter_L_original[ df$splitter_L_original[,"left_child" ]  %in% df$leaf_node_No_L,] 
                                                  } else {
                                                  df$tree #i.e. if save tree here, then N_freq_woN1_nrow=0, already got all parent node levels.
                                                  }
                                                  )              #
List_right <- lapply(List0, 
                     function(df) 
                     if (df$N_freq_woN1_nrow >0) {
                                                  df$splitter_L_original[ df$splitter_L_original[,"right_child"]  %in% df$leaf_node_No_L,]
                                                  } else {
                                                  df$tree  #i.e. if save tree here, then N_freq_woN1_nrow=0, already got all parent node levels.
                                                  }
                                                  )              #  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#  
#step 3.3 add logic sign (subgroup defintion), and position with parent node (TRUE/FALSE),                                                                             #
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#####                                                                                                                   
#if children node on the left, then parent_left_child =TRUE; if children node on the right, then parent_right_child=TRUE,                                              #  
#add "<=" sign, the top 1, 2 node may not be included in lef/right node, thus if (nrow(df) != 0 )                                                                      #
L_sign_left_pre <- list.zip(List_left, N_freq_woN1_nrow, tree)
List_sign_left  <- lapply(L_sign_left_pre,  
                          function(df) 
                          if (df$N_freq_woN1_nrow >0 & nrow(df$List_left)  != 0 ) {
                                                                                  transform(df$List_left, parent_sign = "<=", parent_left_child ='TRUE', parent_right_child ='FALSE') 
                                                                                  } 
                                                                                  )  #
#add ">"  sign, the top 1, 2 node may not be included in lef/right node, thus if (nrow(df) != 0 )                                                                      #
L_sign_right_pre <- list.zip(List_right, N_freq_woN1_nrow, tree)
List_sign_right <- lapply(L_sign_right_pre, 
                          function(df) 
                          if (df$N_freq_woN1_nrow >0 & nrow(df$List_right) != 0 ) {
                                                                                  transform(df$List_right, parent_sign = ">",  parent_left_child ='FALSE',parent_right_child ='TRUE' )   #
                                                                                  } 
                                                                                  )  #
#save both left parent and right parent node info as 2 dataframes in a list                                                                                            #
#in the combo list leaf_L_N1_c, leaf_L_N1_c$leaf_L_N1 as a label for leaf_L_N1_c$leaf_L_N1_f exist or not                                                              #
List1 <- list.zip(List_sign_left, List_sign_right, N_freq_woN1_nrow, tree)                                                                                                               #
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#####  
#step 3.4 #stack (rbind) the left parent and right parent dataframe, if dataframe=NULL, then only the other part will be added                                         #
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------# 
parent <- lapply(List1, 
                 function(df) 
                 if (df$N_freq_woN1_nrow >0) {
                                              rbind(df$List_sign_left, df$List_sign_right)
                                              } else {
                                              df$tree  
                                              }
                                              )                                                                                 #
#already add leading 0 to node number, otherwise when ranking an only be child node rather than parent node,                                                           #
#without leading 0, when ranking, it will be node-1> node-11 > node-12 > ... node-19 > node2                                                                           #
#with    leading 0, when ranking, it will be node-01> node-02 > ...    > ... node-09 > node-10                                                                         #
parent_sort_pre <- list.zip(parent, N_freq_woN1_nrow, tree)      
parent_sort  <- lapply(parent_sort_pre, 
                       function(df)  
                       if (df$N_freq_woN1_nrow >0) {
                                                   df$parent[with(df$parent, order(df$parent[1])), ]
                                                   } else {
                                                   df$tree   
                                                   }
                                                   )                   #order by parent_node_1                                                  #
######################################################################################################################################################################  
######################################################################################################################################################################  
#SUBSECTION 4: [IDENTIFY] parent node # by the current node# of this "parent level"	
######################################################################################################################################################################  
######################################################################################################################################################################  
  if (parent_level==1){
  ######################################################################################################################################################################  
  #Subsection 4.1 lowest parent level of leaves, direct parent!
  ######################################################################################################################################################################  
  #rename the "df" column to the "LEAF" columns to prepare for merging by LEAF for parent level 1.
  #merge LEAF info with parent_sort first! row # are the same, so can merge first!
  List_N0 <- list.zip(parent_sort, LEAF_L_ori ) 
  #To add the leaves for the spliters, the following line to merge is OK, as the original tree dataframe is ordered by number, so is LEAF.
  #deepest nodes will have two rows (leaves) matching to them
  List_N4 <- lapply(List_N0,  function(df)  merge(data.frame(df$parent_sort, row.names=NULL),  data.frame(df$LEAF_L_ori, row.names=NULL), by = 0, all = TRUE)[-1] %>% mutate(Child_node = LEAF ) ) 
} else {
  
  ######################################################################################################################################################################  
  #Subsection 4.3   other parent level of leaves, i.e. [1, max number of leaves]  need to add previous parent node info, 
  ######################################################################################################################################################################  
  #rename the "df" column to the "LEAF" columns to prepare for merging by LEAF for parent level 1.
  leaf_Cnode_FREQ <- lapply(leaf_node_freq_woN1, function(df) {names(df)[names(df) == 'df' ] <- paste('Child_node');  df}) 
  #--------------------------------------------------------------------------------------------------------------------------------------  
  #step 1.  #merge node frequency info with parent_sort since unlike parent_level=1 where parent_sort has the same row # with LEAF info
  #--------------------------------------------------------------------------------------------------------------------------------------  
  List_N0 <- list.zip(parent_sort, leaf_Cnode_FREQ,  N_freq_woN1_nrow, tree)
  #merge Child_node frequency from previous run with current run(parent_sort)!!! 
  #row # are the same, so can merge without worrying about order for now!
  List_N1 <- lapply(List_N0,  
                    function(df)  
                    if (df$N_freq_woN1_nrow >0) {
                                                merge(data.frame(df$parent_sort, row.names=NULL),  data.frame(df$leaf_Cnode_FREQ, row.names=NULL), by = 0, all = TRUE)[-1]
                                                } else {
                                                  df$tree   
                                                }
                                                )   
  #--------------------------------------------------------------------------------------------------------------------------------------  
  #step 2.  subsetting, SORT
  #--------------------------------------------------------------------------------------------------------------------------------------  
  parent_bothC_pre <- list.zip (List_N1, N_freq_woN1_nrow, tree)
      #--------------------------------------------- 
      #step 2.1.  subsetting node w/ Freq==1, 
      #--------------------------------------------- 
      #when level number is high, node w/ Freq==1 may not exist
      parent_bothC_1 <- lapply(parent_bothC_pre, 
                             function(df)   
                             if (df$N_freq_woN1_nrow >0) {
                                                          df$List_N1 [which( (df$List_N1)$Freq==1 ),]
                                                         } else  {
                                                           df$tree   
                                                         }
                                                         )
      #--------------------------------------------- 
      #step 2.2.  subsetting node w/ Freq!=1, 
      #---------------------------------------------  
      #subsetting node w/ Freq!=1, i.e. with nodes with 2 child nodes (LEAVES/subgroups) 
      parent_bothC_gt1 <- lapply(parent_bothC_pre, 
                                 function(df)   
                                 if (df$N_freq_woN1_nrow >0) {
                                                               df$List_N1 [which( (df$List_N1)$Freq!=1 ),] 
                                                             } else  {
                                                               df$tree   
                                                             }
                                                             ) 
      
      #replicate rows the number of times specified in the column "Freq", then correct value in "Freq", which should be 1
      parent_bothC_gt1_D_pre <- list.zip(parent_bothC_gt1, N_freq_woN1_nrow, tree)
      parent_bothC_gt1_D <- lapply(parent_bothC_gt1_D_pre,  
                                   function(df)  
                                   if (df$N_freq_woN1_nrow >0) {
                                                               df$parent_bothC_gt1[rep(row.names(df$parent_bothC_gt1), (df$parent_bothC_gt1)$Freq), 1:ncol(df$parent_bothC_gt1)]  %>% dplyr::mutate(Freq= 1)
                                                               } else  {
                                                                 df$tree   
                                                               }
                                                               )
      #--------------------------------------------- 
      #step 2.3.  Combine nodes Freq!=1 and Freq==1 
      #--------------------------------------------- 
  List_bothC <- list.zip(parent_bothC_1, parent_bothC_gt1_D, N_freq_woN1_nrow, tree)
  #although parent_bothC_1 comes originally from List_N0 which is sorted,  NEED TO RESORT since stacking for nodes freq==1 and !=1 !!!
  parent_bothC <- lapply(List_bothC,  
                         function(df)  
                         if (df$N_freq_woN1_nrow >0) {
                                                      rbind(df$parent_bothC_1, df$parent_bothC_gt1_D)
                                                     } else {
                                                             df$tree
                                                             }
                                                             )
  #SORT by Child_node! key step to merge accurately!!!
  parent_bothC_pre <- list.zip(parent_bothC, N_freq_woN1_nrow, tree)
  parent_bothC_sort <- lapply(parent_bothC_pre, 
                              function(df)  
                              if (df$N_freq_woN1_nrow >0) {
                                                          df$parent_bothC[with(df$parent_bothC, order((df$parent_bothC)$Child_node)), ]
                                                          } else {
                                                                  df$tree   
                                                                  }
                                                                  )  #order by parent_node_1 
  #--------------------------------------------------------------------------------------------------------------------------------------  
  #step 3. generating multiple rows of node-01 as needed
  #-------------------------------------------------------------------------------------------------------------------------------------- 
  #combine preparing for merge
  List_N2 <- list.zip(parent_bothC_sort, leaf_L_N1_co, leaf_node_freq, leaf_node_freq_woN1, N_freq_woN1_nrow, tree)                                                                                                               #
  #--------------------------------------------------------------------------------------------------------------------------------------  
  #step 4 added node-01 as needed, replace NA by node nubmer and Freq 1 for node-01                                          #
  #-------------------------------------------------------------------------------------------------------------------------------------- 
  #add node-01 info
  List_N2_Node1 <- 
    lapply(List_N2,  
          function(df)  
                         if (df$N_freq_woN1_nrow >0 & identical(df$leaf_node_freq, df$leaf_node_freq_woN1)==TRUE)
                  {
                   df$parent_bothC_sort
                  } else if (df$N_freq_woN1_nrow >0 & identical(df$leaf_node_freq, df$leaf_node_freq_woN1)==FALSE)
                  {
                   dplyr::bind_rows(df$leaf_L_N1_co, df$parent_bothC_sort) 
                  } else {
                    df$tree 
                  }
                  )       
  #keep node-01 info when necessary, i.e. parent_sign, parent_left_child, or parent_right_child != NA, drop if =NA
  List_N2_Node1_T_pre <- list.zip(List_N2_Node1, N_freq_woN1_nrow, tree)
  
  List_N2_Node1_T  <- lapply (List_N2_Node1_T_pre, 
                             function(df)  
                             if (df$N_freq_woN1_nrow >0) {
                                                         df$List_N2_Node1 <- df$List_N2_Node1[!(is.na ((df$List_N2_Node1)$parent_sign)==TRUE ),] 
                                                         } else {
                                                                 df$tree  
                                                                 }
                                                                 )
  
  #connect original LEAF node, and parent_level info, prepare for merge, then replace "NA": of Freq column by 1
  List_N3 <- list.zip(List_N2_Node1_T, LEAF_L_ori,  N_freq_woN1_nrow, tree)
  List_N4 <- 
    lapply(List_N3 ,  
           function(df) 
           if (df$N_freq_woN1_nrow >0) {
                                       merge(data.frame(df$List_N2_Node1_T, row.names=NULL),  data.frame(df$LEAF_L_ori, row.names=NULL), by = 0 , all = TRUE)[-1] %>% dplyr::mutate(Freq = tidyr::replace_na(Freq, 1))
                                       } else {
                                              df$tree 
                                              }
                                              ) 

}
########################################################################################################################################################################  
########################################################################################################################################################################  
#SUBSECTION 5: [OUTPUT] the parent level info 
########################################################################################################################################################################  
########################################################################################################################################################################  
#
parent_condition_pre <- list.zip (List_N4, N_freq_woN1_nrow, tree)
parent_condition <- lapply(parent_condition_pre, 
                           function(df) 
                           if (df$N_freq_woN1_nrow >0) {
                                                        transform(df$List_N4, condition = paste((df$List_N4)$split_variable, (df$List_N4)$parent_sign, (df$List_N4)$split_value, sep="") )
                                                        } else {
                                                          df$tree 
                                                        }
                                                        ) 


#rename to show parent_level to preparing merging with next parent level                                                                                                             #
parent_rename_pre <- list.zip (parent_condition, N_freq_woN1_nrow, tree)
parent_rename <- lapply(parent_rename_pre, 
                        function(df) 
                        if (df$N_freq_woN1_nrow >0) {
                                                     {names(df$parent_condition)[names(df$parent_condition) == 'parent_left_child' ] <- paste('parent_left_child',  parent_level, sep="_");            #   
                                                      names(df$parent_condition)[names(df$parent_condition) == 'parent_right_child'] <- paste('parent_right_child', parent_level, sep="_");            # 
                                                      names(df$parent_condition)[names(df$parent_condition) == 'parent_sign']        <- paste('parent_sign',        parent_level, sep="_");            #             
                                                      names(df$parent_condition)[names(df$parent_condition) == 'split_variable']     <- paste('parent_split_var',   parent_level, sep="_");            # 
                                                      names(df$parent_condition)[names(df$parent_condition) == 'split_value']        <- paste('parent_split_val',   parent_level, sep="_");            # 
                                                      names(df$parent_condition)[names(df$parent_condition) == 'node']               <- paste('node',               parent_level, sep="_");            #
                                                      names(df$parent_condition)[names(df$parent_condition) == 'Child_node']         <- paste('Child_node',         parent_level, sep="_");            #
                                                      names(df$parent_condition)[names(df$parent_condition) == 'Freq']               <- paste('Freq',               parent_level, sep="_");            #
                                                      names(df$parent_condition)[names(df$parent_condition) == 'left_child']         <- paste('left_child',         parent_level, sep="_");            #
                                                      names(df$parent_condition)[names(df$parent_condition) == 'right_child']        <- paste('right_child',        parent_level, sep="_"); 
                                                      names(df$parent_condition)[names(df$parent_condition) == 'condition']          <- paste('condition',          parent_level, sep="_"); 
                                                      df$parent_condition}
                                                      } else {
                                                              df$tree
                                                             }
                                                             )                                                                                                                        #  


    
#remove non related columns                                                                                                                                            #
parent_final_pre <- list.zip(parent_rename, N_freq_woN1_nrow, tree)
parent_final <- lapply(parent_final_pre,
                       function(df) 
                       if (df$N_freq_woN1_nrow >0) {
                                                    df$parent_rename[, -which(names(df$parent_rename) %in% c("is_leaf","samples", "avg_Y", "avg_W", "k", #iCF
                                                                                                                                                    "b", #oneCF
                                                                                                                                                    "HTE_P_cf"))]
                                                   } else {
                                                           df$tree                       
                                                          }
                                                          )                      #
########################################################################################################################################################################  
# Step 5.1 output final result depending on parent level, combine with previous output when appropriate                                                                #
########################################################################################################################################################################  
if (parent_level==1){                                                                                                                                                  #
                    parent_final <- parent_final                                                                                                                       #
                    } else {                                                                                                                                           #
                            #------------------------------------------------------------------------------------------------------------------------------------------#  
                            # combine with previoous output  (i.e. current tree) if parent_level > 1 & non-node-01 still exists (N_freq_woN1_nrow >0)                  #
                            #------------------------------------------------------------------------------------------------------------------------------------------#  
                            List_N5 <- list.zip(parent_final, N_freq_woN1_nrow, tree)                                                                                                   
                            parent_final <- lapply(List_N5,  
                                                   function(df) 
                                                   if (df$N_freq_woN1_nrow >0) {
                                                                               left_join(df$parent_final, df$tree, by = "LEAF" )
                                                                               } else {
                                                                                       df$tree  
                                                                                       }
                                                                                       ) #12/17/2020 updates!!!                                         
                          }                                                                                                                                                                      
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#  


return(parent_final)


  

  }


