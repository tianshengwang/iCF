
#----------------------------------------------------------------------------
#7/29/2023 blocked the following, seems useless, will delete later.
#----------------------------------------------------------------------------
#HTE_P_extract_M <- function(iCF_D_M, tree_depth){
#  HTE_P_cf_M <- iCF_D_M$HTE_P_cf_M
#  HTE_P_cf_M$tree_depth <- tree_depth
#  return(HTE_P_cf_M)
#}

#SF_extract_M <- function(iCF_D_M, tree_depth){
#  SF_M <- iCF_D_M$SF_ID_LIST
#  return(SF_M)
#}

#combine all df function:
#COMBO_SAVE <- function(List,  Listname){
#  allonelist = do.call(rbind, List)
#record allbesttree data
#  write.table(allonelist,file=paste(Listname, depth, treeNo, iterationNo, Ntrain, tree_depth, intTRUE, ".csv", sep="_"),row.names=FALSE,col.names=T,sep=",",append=TRUE)
#save the list of frames
#  save(List,             file=paste(Listname, depth, treeNo, iterationNo, Ntrain, tree_depth, intTRUE, ".Rda", sep="_"))
#}

#FOREST_LIST <- function(dataset) {



#  vars_forest =c("X1","X2","X3","X4","X5","X6","X8","X9","X10") 
#  vars_IV=c("X7")
#  X<-dataset[,vars_forest]
#  Y<-dataset[,"Y"]
#  W<-dataset[,"a"]
#  Z<-dataset[,vars_IV]  
#  forest.list <- list(X = X, Y = Y, W = W, Z = Z)
#  return(forest.list)	    
#}
#------------------------------------------------------------------
#7/29/2023, not sure the following is useful, may delete later
#------------------------------------------------------------------
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

#' This function Impute split values for splitters of D3 voted tree from D4 synthetic tree if top 3 nodes are identical
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



#' FL_Model_Selection
#' 
#' Function that does backward selection of best model, works for both continuous (partial F-test) and binary outcome (likelihood ratio test)
#' @param P_basicvsg2 P-value for partial F-test or LRT for main effect model vs D2 interaction term model
#' @param P_g2vsg23   P-value for partial F-test or LRT for D2 interaction term model vs D3 and D2 interaction terms model 
#' @param P_g23vsg234 P-value for partial F-test or LRT for D3 and D2 interaction term model vs D4, D3 and D2 interaction terms model  
#' @param siglevel predefined significance level for P-value
#'  
#' @return the data for group lasso
#' 
#' @export
FL_Model_Selection <- function(P_basicvsg2, P_g2vsg23, P_g23vsg234, siglevel){
  if ( 
    #scenario 1: all significant (keep G4)
    max(P_basicvsg2, P_g2vsg23, P_g23vsg234 ) < siglevel ) {
    grpFanova_pick_D = "W:G4x"
  } else if ( 
    #scenario 2: G23 vs G4 not significant (drop G4), remaining 2 both significant (keep G3) 
    max(P_basicvsg2, P_g2vsg23) < siglevel & P_g23vsg234 >= siglevel ) {
    grpFanova_pick_D = "W:G3x"
  } else if (
    #scenario 3: neither G23 vs G4 not nor G2 vs G23 are significant (dropped G3), main vs G2 significant (keep G2)
    ( P_basicvsg2 < siglevel & min(P_g2vsg23, P_g23vsg234) >= siglevel |
      #scenario 4: both G23 vs G4 and G2 vs basic are significant but G2 vs G23 not significant (dropped G3)
      P_basicvsg2 < siglevel & P_g2vsg23 >= siglevel &  P_g23vsg234 < siglevel ) ) {
    grpFanova_pick_D = "W:G2x"
  } else if ( 
    #scenario 5: main vs G2 not significant
    P_basicvsg2 >= siglevel  ) { 
    grpFanova_pick_D = "NA"
  }
  return(grpFanova_pick_D)
}



#' FINALDECI, blocked on 7/29/2023
#' 
#' Function that does shows the final subgroup decisions 
#' @param grp_pick 
#' @param V_D4_subgroup   
#'  
#' @return the final subgroup decision,
#' 
#' @export
#FINALDECI<- function(grp_pick, V_D4_subgroup,  V_D3_subgroup, V_D2_subgroup){
#  CFmethod_prefix  = "vote_D"
#  CFmethod_postfix = "_subgroup"
#  V_list=list(D2= V_D2_subgroup, D3= V_D3_subgroup,D4= V_D4_subgroup)
#  if (grp_pick != "NA"){
#    Deci_Final        <-  V_list[[as.numeric( stringr::str_sub(grp_pick,     4,4)) - 1]]$majority
#  } else {
#    Deci_Final = "NA"
#  }
#  return(Deci_Final)
#}




#7/29/2023 removed the following, useless
#------------------------------------------------------------------------------------------------------------
# IV. ACCURACY (=TRUE OR NOT) OF EACH TYPE OF TREE AND FINAL SUBGROUP DECISION BY DIFFERENT THRESHOLD VALUE 
#------------------------------------------------------------------------------------------------------------
#TEST_DECISION <- function(test_df, true_df){
#  if (identical(test_df, true_df)){ 
#    test.score <- 1
#  } else {
#    test.score <- 0
#  }
#  return(test.score)
#}

#------------------------------------------------------------------------------------------------------------
# IV. Which Depth iCF decision comes from 
#------------------------------------------------------------------------------------------------------------
#WHICH_DECISION <- function(depth_df, final_df){
#  if (identical(depth_df, final_df)){ 
#    test.score <- 1
#  } else {
#    test.score <- 0
#  }
#  return(test.score)
#}

#'function that check tree structure/subgroup are true or not 
#' @param v_D4_tree voted D4 tree
#' @param v_D4_tree_R voted D4 tree without leaves
#' @param Desc_v_D4 decision from D4 tree
#' @param v_D3_tree voted D3 tree
#' @param v_D3_tree_R voted D3 tree without leaves
#' @param Desc_v_D3 decision from D3 tree
#' @param v_D2_tree voted D2 tree
#' @param v_D2_tree_R voted D2 tree without leaves
#' @param Desc_v_D2 decision from D2 tree 
#' @return subgroup in list format
#' 
#' @export
#TorF <- function(v_D4_tree, v_D4_tree_R, Desc_v_D4, v_D3_tree, v_D3_tree_R, Desc_v_D3, v_D2_tree, v_D2_tree_R, Desc_v_D2, truth.list){
#  TF_D4_T     <- ifelse (identical(v_D4_tree$majority,      truth.list$tree_true), 1, 0)
#  TF_D4_T_N1  <- ifelse (identical(v_D4_tree$majority[1,],  truth.list$tree_true_N1), 1, 0)   #get from tree list, cannot check from subgroup or tree_r list
#  TF_D4_T_N123<- ifelse (identical(v_D4_tree$majority[1:3,],truth.list$tree_true_N123), 1, 0) #get from tree list, cannot check from subgroup or tree_r list
#  TF_D4_T_r   <- ifelse (identical(v_D4_tree_R$majority,    truth.list$tree_true_r), 1, 0)
#  TF_D4_SG    <- ifelse (identical(Desc_v_D4,               truth.list$tree_true_subgroup), 1, 0)
#  
#  TF_D3_T     <- ifelse (identical(v_D3_tree$majority,      truth.list$tree_true), 1, 0)
#  TF_D3_T_N1  <- ifelse (identical(v_D3_tree$majority[1,],  truth.list$tree_true_N1), 1, 0)   #get from tree list, cannot check from subgroup or tree_r list
#  TF_D3_T_N123<- ifelse (identical(v_D3_tree$majority[1:3,],truth.list$tree_true_N123), 1, 0) #get from tree list, cannot check from subgroup or tree_r list
#  TF_D3_T_r   <- ifelse (identical(v_D3_tree_R$majority,    truth.list$tree_true_r), 1, 0)
#  TF_D3_SG    <- ifelse (identical(Desc_v_D3,               truth.list$tree_true_subgroup), 1, 0)
#  
#  TF_D2_T     <- ifelse (identical(v_D2_tree$majority,      truth.list$tree_true), 1, 0)
#  TF_D2_T_N1  <- ifelse (identical(v_D2_tree$majority[1,],  truth.list$tree_true_N1), 1, 0)   #get from tree list, cannot check from subgroup or tree_r list
#  TF_D2_T_N123<- ifelse (identical(v_D2_tree$majority[1:3,],truth.list$tree_true_N123), 1, 0) #always wrong (i.e 0) when truth is "deeper tree"
#  TF_D2_T_r   <- ifelse (identical(v_D2_tree_R$majority,    truth.list$tree_true_r), 1, 0)
#  TF_D2_SG    <- ifelse (identical(Desc_v_D2,               truth.list$tree_true_subgroup), 1, 0)
#  
#  return(
#    list(  TF_D4_T     = TF_D4_T, 
#           TF_D4_T_N1  = TF_D4_T_N1,
#           TF_D4_T_N123= TF_D4_T_N123,
#           TF_D4_T_r   = TF_D4_T_r ,
#           TF_D4_SG    = TF_D4_SG ,
#           
#           TF_D3_T     = TF_D3_T,
#           TF_D3_T_N1  = TF_D3_T_N1,
#           TF_D3_T_N123= TF_D3_T_N123, 
#           TF_D3_T_r   = TF_D3_T_r ,
#           TF_D3_SG    = TF_D3_SG,
#           
#           TF_D2_T     = TF_D2_T, 
#           TF_D2_T_N1  = TF_D2_T_N1 ,
#           TF_D2_T_N123= TF_D2_T_N123 , 
#           TF_D2_T_r   = TF_D2_T_r ,
#           TF_D2_SG    = TF_D2_SG  
#    )
#  )
#}


#DISCOVERY <- function(N_method, N_truth, TF){
#  discovery="NA"
#  if (N_method < N_truth){
#    discovery="typeII"
#  } else if (N_method > N_truth ) {
#    discovery="typeI"
#  } else if (N_method == N_truth & TF==1) {
#    discovery="accurate"
#  }
#  return(discovery)
#}



#' Function that check if root node of tree structure is accurate  
#' @param decision the final decision obtained from the decision path
#' @param vote_D4_sg the list of D4 vote subgroup results
#' @param vote_D3_sg the list of D3 vote subgroup results
#' @param vote_D2_sg the list of D2 vote subgroup results
#' @param vote_D4_t the list of D4 vote tree results
#' @param vote_D3_t the list of D3 vote tree results
#' @param vote_D2_t the list of D2 vote tree results
#' @param tree_true the truth of tree structure 
#' 
#' @return if root node of tree structure is accurate   
#' 
#' @export

#DECI_treeN1_mCF_TF <- function(decision, vote_D4_sg, vote_D3_sg, vote_D2_sg, vote_D4_t, vote_D3_t, vote_D2_t, tree_true){
#  #first identify the decision comes from which depth tree 
#  if (identical( decision, vote_D4_sg$majority ) ) {
#    DECI_treeN1 <<- vote_D4_t$majority[1, "split_variable"] 
#  } else if (identical( decision, vote_D3_sg$majority ) ) {
#    DECI_treeN1 <<- vote_D3_t$majority[1, "split_variable"] 
#  } else if (identical( decision, vote_D2_sg$majority ) ) {
#    DECI_treeN1 <<- vote_D2_t$majority[1, "split_variable"] 
#  } else if ( identical( decision, "NA" ) ){
#    DECI_treeN1 <<- "NA"
#  }
  #second check if top node is accurate
#  if (identical(tree_true, "NA")==F) {
#    TorF= ifelse (identical(DECI_treeN1, tree_true[1, "split_variable"]), 1, 0)
#  } else if (tree_true == "NA") {
#    TorF= ifelse (identical(DECI_treeN1, tree_true), 1, 0)
#  }
#  return(TorF)
#}

#' Function that check if the tree of the subgroup decision has the correct top node 
#' @param decision the subgroup decision 
#' @param vote_sg 
#' @param vote_t
#' @param tree_true true tree structure
#' 
#' @return subgroup decision   
#' 
#' @export
#DECI_treeN1_TREE_TF <- function(decision, vote_sg, vote_t, tree_true){
#  if ( identical(decision, "NA")==F ) {
#    DECI_treeN1 <<- vote_t[1, "split_variable"] 
#  } else if (identical( decision, "NA" ) ){
#    DECI_treeN1 <<- "NA"
# }
#  if ( identical(tree_true, "NA")==F) {
#    TorF= ifelse (identical(DECI_treeN1, tree_true[1, "split_variable"]), 1, 0)
#  } else if (tree_true == "NA") {
#    TorF= ifelse (identical(DECI_treeN1, tree_true), 1, 0)
#  }
#  return( TorF )
#}




#' Function that run cross-validation for group lasso (grplasso) using a formula (formula_grplasso) designed to CF, 
#' BE CAREFUL!!! Need to revise for classification.
#' @param dat the dataset 
#' @param lambda_grplasso a grid of penalty parameter from 0 to lambda_max
#' 
#' @return the list of cross-validation results and the lambda leading to minimor. 
#' 
#' @export

cv.grplasso <- function(dat, lambda_grplasso, formula_grplasso, contrast){
  
  env=global_env() #Advanced R, P477!!!34BF98CB|
  gl_cv_results <- list()
  
  
  for (j in 1:length(lambda_grplasso)){
    
    
    mse_per_fold      <- list()
    accuracy_per_fold <- list()
    
    for(i in 1:5){
      cv_indices <- caret::createFolds(y=dat$Y, k=5)
      train_data <- dat[-cv_indices[[i]],]
      test_data  <- dat[cv_indices[[i]],]
      
      
      if (length(unique(dat$Y)) >=8){
        #==================================================
        # Create linear regression model for Group Lasso
        #==================================================
        grplasso_fit <-  grplasso::grplasso(formula = formula_grplasso, data =  train_data, contrasts = contrast,  center = TRUE, standardize = TRUE, lambda = lambda_grplasso[j],
                                            model = grplasso::LinReg()
        )
        pred <- predict(grplasso_fit, newdata = test_data)
        test_data$predict_gl <- predict(grplasso_fit, newdata=test_data)
        mse_per_fold[[i]] <- caret::postResample(pred = test_data$predict_gl, obs = test_data$Y)
      } else if (length(unique(dat$Y)) == 2 ){
        #==================================================
        # Create logistic regression model for Group Lasso
        #==================================================
        grplasso_fit <-  grplasso::grplasso(formula = formula_grplasso, data =  train_data, contrasts = contrast,  center = TRUE, standardize = TRUE, lambda = lambda_grplasso[j],
                                            model = grplasso::LogReg()
        ) 
        test_data$predict_gl <- as.numeric(predict(grplasso_fit, newdata=test_data, type="response"))
        hist(test_data$predict_gl)
        test_data <- test_data %>% dplyr:: mutate(predict_Y = factor(ifelse(predict_gl>0.5, 1, 0) ),
                                                  Y=factor(Y))
        
        table(test_data$predict_Y)
        
        ConfusionM <- caret::confusionMatrix(data = test_data$predict_Y,
                                             reference = test_data$Y,
                                             positive = "0" #Medicare DPP4i vs SU; 0=event, 1=no event so that larger outcome is better
        )
        accuracy_per_fold[[i]] <- cbind("lambda"=lambda_grplasso[j], 
                                        "Accuracy"=ConfusionM[[3]][1])
      }
      
    }
    ##############################  
    # i iteration ends
    ##############################  
    
    
    if (length(unique(dat$Y)) >=8){
      # Bind together, add MSE
      mse_all_folds <- as.data.frame(do.call("rbind", mse_per_fold))
      mse_all_folds <- cbind(mse_all_folds, "MSE"=(mse_all_folds[,"RMSE"])^2)
      # Get mean and SE MSE
      gl_cv_results[[i]] <- data.frame("j"=j,
                                       "lambda"=lambda_grplasso[j],
                                       "CV_MSE"=mean(mse_all_folds[,"MSE"]),
                                       "CV_MSE_SE"=sd(mse_all_folds[,"MSE"]))
    } else if (length(unique(dat$Y)) ==2){
      # Bind together, add Accuracy
      accuracy_all_folds <- as.data.frame(do.call("rbind", accuracy_per_fold))
      # Get mean accuracy
      gl_cv_results[[j]] <- data.frame("j"= j , 
                                       "lambda"=lambda_grplasso[j], 
                                       CV_Accuracy=mean(accuracy_all_folds$Accuracy)
      )
    }
    
  } 
  ##############################  
  # j iteration ends
  ##############################
  #combine all results in the list into a dataframe
  gl_cv_results <- do.call("rbind", gl_cv_results)
  
  
  if (length(unique(dat$Y)) >=8){
    #continouous outcome: the lambda leads to minimum MSE
    lambda.min = lambda_grplasso[ which(gl_cv_results$CV_MSE == min(gl_cv_results$CV_MSE)) ]
  } else if  (length(unique(dat$Y)) ==2 ) {
    #binary outcome: the lambda leads to maximum accuracy
    lambda.min = lambda_grplasso[ which(gl_cv_results$CV_Accuracy == max(gl_cv_results$CV_Accuracy)) ]
  }
  
  
  return(list(cv.grplasso = gl_cv_results, lambda.min = lambda.min))
}



#' DELTA_Y_DATA_mCF
#' 
#' Function that extract ATE, ATT, and crude Delta Y and the coeffecient of GLM for ggploting error bar of Delta Y for both training and testing data. 
#' Run SUBSETTING, ANA_SUBGROUP, and EXTRACT_SUB_ANA 3 functions
#' Forest could be oneCF, iCF, or iCFv
#' @param Deci_D4 the voted decision from D4 forest 
#' @param Deci_D3 the voted decision from D3 forest 
#' @param Deci_D2 the voted decision from D2 forest
#' @train Training dataset
#' @test Testing dataset
#' @methodL method label, "oneCF", "iCF", or "iCFv"
#' @return The list of extracted ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) for both training and testing data, and model df of coeffecient for treatment for ATE, ATT, and crude 
#' 
#' @export

DELTA_Y_DATA_mCF <- function(Deci_D4, Deci_D3, Deci_D2, traindata, testdata, methodL){
  #make sure Fx work in another Fx: call global environment
  env=global_env() 
  
  SG_D4_tr <- SUBSETTING(Deci_D4, traindata); SG_D3_tr <- SUBSETTING(Deci_D3, traindata); SG_D2_tr <- SUBSETTING(Deci_D2, traindata)
  SG_D4_te <- SUBSETTING(Deci_D4, testdata);  SG_D3_te <- SUBSETTING(Deci_D3, testdata);  SG_D2_te <- SUBSETTING(Deci_D2, testdata)
  #checking if sum of each subgroup # match total N
  nrow( bind_rows(SG_D4_tr, .id = "column_label") ); nrow( bind_rows(SG_D3_tr, .id = "column_label") ); nrow( bind_rows(SG_D2_tr, .id = "column_label") )
  nrow( bind_rows(SG_D4_te, .id = "column_label") ); nrow( bind_rows(SG_D3_te, .id = "column_label") ); nrow( bind_rows(SG_D2_te, .id = "column_label") )
  
  SG_D4_y_tr <- lapply(SG_D4_tr, ANA_SUBGROUP); SG_D3_y_tr <- lapply(SG_D3_tr, ANA_SUBGROUP); SG_D2_y_tr <- lapply(SG_D2_tr, ANA_SUBGROUP)
  SG_D4_y_te <- lapply(SG_D4_te, ANA_SUBGROUP); SG_D3_y_te <- lapply(SG_D3_te, ANA_SUBGROUP); SG_D2_y_te <- lapply(SG_D2_te, ANA_SUBGROUP)
  
  #----------------------------
  # Extract Delta Y and CI
  #----------------------------
  #subsetting nested list
  sg_D4_y_tr <- lapply(SG_D4_y_tr, "[", 1);                        sg_D3_y_tr <- lapply(SG_D3_y_tr, "[", 1);                        sg_D2_y_tr <- lapply(SG_D2_y_tr, "[", 1)	
  sg_D4_y_te <- lapply(SG_D4_y_te, "[", 1);                        sg_D3_y_te <- lapply(SG_D3_y_te, "[", 1);                        sg_D2_y_te <- lapply(SG_D2_y_te, "[", 1)	
  #"unnested" the dataframe which was in a list nested in a list
  sg_D4_y_tr_un <- lapply(sg_D4_y_tr, function(df) df$subgroup);   sg_D3_y_tr_un <- lapply(sg_D3_y_tr, function(df) df$subgroup);   sg_D2_y_tr_un <- lapply(sg_D2_y_tr, function(df) df$subgroup)
  sg_D4_y_te_un <- lapply(sg_D4_y_te, function(df) df$subgroup);   sg_D3_y_te_un <- lapply(sg_D3_y_te, function(df) df$subgroup);   sg_D2_y_te_un <- lapply(sg_D2_y_te, function(df) df$subgroup)
  #extract Delta Y                           
  Dlt_Y_D4_tr<-EXTRACT_SUB_ANA(sg_D4_y_tr_un,4,"train", methodL);  Dlt_Y_D3_tr<-EXTRACT_SUB_ANA(sg_D3_y_tr_un,3,"train", methodL);  Dlt_Y_D2_tr<-EXTRACT_SUB_ANA(sg_D2_y_tr_un,2,"train",methodL)
  Dlt_Y_D4_te<-EXTRACT_SUB_ANA(sg_D4_y_te_un,4,"test",  methodL);  Dlt_Y_D3_te<-EXTRACT_SUB_ANA(sg_D3_y_te_un,3,"test",  methodL);  Dlt_Y_D2_te<-EXTRACT_SUB_ANA(sg_D2_y_te_un,2,"test", methodL)
  #make a dataframe: use rbind.data.frame to avoid rbind()/cbind() conversion from numeric to factor!!!
  Dlt_Y_tr <-  rbind.data.frame(Dlt_Y_D4_tr, Dlt_Y_D3_tr, Dlt_Y_D2_tr) 
  Dlt_Y_te <-  rbind.data.frame(Dlt_Y_D4_te, Dlt_Y_D3_te, Dlt_Y_D2_te) 
  #keep the original order (by SubID, SubDef, Depth)
  Dlt_Y <- plyr::join(Dlt_Y_tr, Dlt_Y_te)#	%>% mutate_at( vars(-SubID, -SubDef, -Depth), as.factor)
  
  #----------------------------
  # Extract Delta Y and CI
  #----------------------------
  #subsetting nested list for 3 models
  sg_D4_y_tr_ate <- lapply(SG_D4_y_tr, "[", 2);  sg_D3_y_tr_ate <- lapply(SG_D3_y_tr, "[", 2);  sg_D2_y_tr_ate <- lapply(SG_D2_y_tr, "[", 2);
  sg_D4_y_te_ate <- lapply(SG_D4_y_te, "[", 2);  sg_D3_y_te_ate <- lapply(SG_D3_y_te, "[", 2);  sg_D2_y_te_ate <- lapply(SG_D2_y_te, "[", 2);
  
  sg_D4_y_tr_att <- lapply(SG_D4_y_tr, "[", 3);  sg_D3_y_tr_att <- lapply(SG_D3_y_tr, "[", 3);  sg_D2_y_tr_att <- lapply(SG_D2_y_tr, "[", 3);
  sg_D4_y_te_att <- lapply(SG_D4_y_te, "[", 3);  sg_D3_y_te_att <- lapply(SG_D3_y_te, "[", 3);  sg_D2_y_te_att <- lapply(SG_D2_y_te, "[", 3);
  
  sg_D4_y_tr_cru <- lapply(SG_D4_y_tr, "[", 4);  sg_D3_y_tr_cru <- lapply(SG_D3_y_tr, "[", 4);  sg_D2_y_tr_cru <- lapply(SG_D2_y_tr, "[", 4);
  sg_D4_y_te_cru <- lapply(SG_D4_y_te, "[", 4);  sg_D3_y_te_cru <- lapply(SG_D3_y_te, "[", 4);  sg_D2_y_te_cru <- lapply(SG_D2_y_te, "[", 4);
  
  #"unnested" the dataframe which was in a list nested in a list for 3 models
  sg_D4_y_tr_ate_un <- lapply(sg_D4_y_tr_ate, function(df) df$model_ate); sg_D3_y_tr_ate_un <- lapply(sg_D3_y_tr_ate, function(df) df$model_ate); sg_D2_y_tr_ate_un <- lapply(sg_D2_y_tr_ate, function(df) df$model_ate); 
  sg_D4_y_te_ate_un <- lapply(sg_D4_y_te_ate, function(df) df$model_ate); sg_D3_y_te_ate_un <- lapply(sg_D3_y_te_ate, function(df) df$model_ate); sg_D2_y_te_ate_un <- lapply(sg_D2_y_te_ate, function(df) df$model_ate); 
  
  sg_D4_y_tr_att_un <- lapply(sg_D4_y_tr_att, function(df) df$model_att); sg_D3_y_tr_att_un <- lapply(sg_D3_y_tr_att, function(df) df$model_att); sg_D2_y_tr_att_un <- lapply(sg_D2_y_tr_att, function(df) df$model_att); 
  sg_D4_y_te_att_un <- lapply(sg_D4_y_te_att, function(df) df$model_att); sg_D3_y_te_att_un <- lapply(sg_D3_y_te_att, function(df) df$model_att); sg_D2_y_te_att_un <- lapply(sg_D2_y_te_att, function(df) df$model_att); 
  
  sg_D4_y_tr_cru_un <- lapply(sg_D4_y_tr_cru, function(df) df$model_cru); sg_D3_y_tr_cru_un <- lapply(sg_D3_y_tr_cru, function(df) df$model_cru); sg_D2_y_tr_cru_un <- lapply(sg_D2_y_tr_cru, function(df) df$model_cru); 
  sg_D4_y_te_cru_un <- lapply(sg_D4_y_te_cru, function(df) df$model_cru); sg_D3_y_te_cru_un <- lapply(sg_D3_y_te_cru, function(df) df$model_cru); sg_D2_y_te_cru_un <- lapply(sg_D2_y_te_cru, function(df) df$model_cru); 
  
  #combine glm model from each subgroup into 1 df
  D4_tr_ate <-bind_rows(sg_D4_y_tr_ate_un, .id = "column_label");   D3_tr_ate <-bind_rows(sg_D3_y_tr_ate_un, .id = "column_label");   D2_tr_ate <-bind_rows(sg_D2_y_tr_ate_un, .id = "column_label");
  D4_te_ate <-bind_rows(sg_D4_y_te_ate_un, .id = "column_label");   D3_te_ate <-bind_rows(sg_D3_y_te_ate_un, .id = "column_label");   D2_te_ate <-bind_rows(sg_D2_y_te_ate_un, .id = "column_label"); 
  
  D4_tr_att <-bind_rows(sg_D4_y_tr_att_un, .id = "column_label");   D3_tr_att <-bind_rows(sg_D3_y_tr_att_un, .id = "column_label");   D2_tr_att <-bind_rows(sg_D2_y_tr_att_un, .id = "column_label");
  D4_te_att <-bind_rows(sg_D4_y_te_att_un, .id = "column_label");   D3_te_att <-bind_rows(sg_D3_y_te_att_un, .id = "column_label");   D2_te_att <-bind_rows(sg_D2_y_te_att_un, .id = "column_label"); 
  
  D4_tr_cru <-bind_rows(sg_D4_y_tr_cru_un, .id = "column_label");   D3_tr_cru <-bind_rows(sg_D3_y_tr_cru_un, .id = "column_label");   D2_tr_cru <-bind_rows(sg_D2_y_tr_cru_un, .id = "column_label");
  D4_te_cru <-bind_rows(sg_D4_y_te_cru_un, .id = "column_label");   D3_te_cru <-bind_rows(sg_D3_y_te_cru_un, .id = "column_label");   D2_te_cru <-bind_rows(sg_D2_y_te_cru_un, .id = "column_label"); 
  
  all_list <- list(Dlt_Y=Dlt_Y,
                   D4_tr_ate = D4_tr_ate,  D3_tr_ate = D3_tr_ate,   D2_tr_ate = D2_tr_ate,
                   D4_te_ate = D4_te_ate,  D3_te_ate = D3_te_ate,   D2_te_ate = D2_te_ate,
                   
                   D4_tr_att = D4_tr_att,  D3_tr_att = D3_tr_att,   D2_tr_att = D2_tr_att,
                   D4_te_att = D4_te_att,  D3_te_att = D3_te_att,   D2_te_att = D2_te_att,
                   
                   D4_tr_cru = D4_tr_cru,  D3_tr_cru = D3_tr_cru,   D2_tr_cru = D2_tr_cru,
                   D4_te_cru = D4_te_cru,  D3_te_cru = D3_te_cru,   D2_te_cru = D2_te_cru)
  
  return(all_list)
}




#' Function that extract ATE, ATT, and crude Delta Y and the coeffecient of GLM for ggploting error bar of Delta Y for both training and testing data. 
#' Run SUBSETTING, ANA_SUBGROUP, and EXTRACT_SUB_ANA 3 functions
#' for Interaction Tree  
#' @param Deci the decision
#' @param traindata Training dataset
#' @param testdata Testing dataset
#' @param methodL method label, "IT"
#' 
#' @return The list of extracted ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) for both training and testing data, and model df of coeffecient for treatment for ATE, ATT, and crude 
#' 
#' @export
#' 

DELTA_Y_DATA <- function(Deci, traindata, testdata, methodL){
  #make sure Fx work in another Fx: call global environment
  env=global_env() 
  
  SG_D_tr <- SUBSETTING(Deci, traindata); 
  SG_D_te <- SUBSETTING(Deci, testdata);  
  #checking if sum of each subgroup # match total N
  nrow( bind_rows(SG_D_tr, .id = "column_label") );
  nrow( bind_rows(SG_D_te, .id = "column_label") );
  
  SG_D_y_tr <- lapply(SG_D_tr, ANA_SUBGROUP);
  SG_D_y_te <- lapply(SG_D_te, ANA_SUBGROUP); 
  
  #----------------------------
  # Extract Delta Y and CI
  #----------------------------
  #subsetting nested list
  sg_D_y_tr <- lapply(SG_D_y_tr, "[", 1);                       
  sg_D_y_te <- lapply(SG_D_y_te, "[", 1);                      
  #"unnested" the dataframe which was in a list nested in a list
  sg_D_y_tr_un <- lapply(sg_D_y_tr, function(df) df$subgroup);   
  sg_D_y_te_un <- lapply(sg_D_y_te, function(df) df$subgroup);  
  #extract Delta Y, Depth=0 ,i.e. not multiple depth!!!                           
  Dlt_Y_D_tr<-EXTRACT_SUB_ANA(sg_D_y_tr_un, 0,"train", methodL); 
  Dlt_Y_D_te<-EXTRACT_SUB_ANA(sg_D_y_te_un, 0,"test",  methodL);  
  #make a dataframe: use rbind.data.frame to avoid rbind()/cbind() conversion from numeric to factor!!!
  Dlt_Y_tr <-  rbind.data.frame(Dlt_Y_D_tr) 
  Dlt_Y_te <-  rbind.data.frame(Dlt_Y_D_te) 
  #keep the original order (by SubID, SubDef, Depth)
  Dlt_Y <- plyr::join(Dlt_Y_tr, Dlt_Y_te)#	%>% mutate_at( vars(-SubID, -SubDef, -Depth), as.factor)
  
  #----------------------------
  # Extract Delta Y and CI
  #----------------------------
  #subsetting nested list for 3 models
  sg_D_y_tr_ate <- lapply(SG_D_y_tr, "[", 2); 
  sg_D_y_te_ate <- lapply(SG_D_y_te, "[", 2);  
  
  sg_D_y_tr_att <- lapply(SG_D_y_tr, "[", 3);  
  sg_D_y_te_att <- lapply(SG_D_y_te, "[", 3);  
  
  sg_D_y_tr_cru <- lapply(SG_D_y_tr, "[", 4); 
  sg_D_y_te_cru <- lapply(SG_D_y_te, "[", 4);  
  
  #"unnested" the dataframe which was in a list nested in a list for 3 models
  sg_D_y_tr_ate_un <- lapply(sg_D_y_tr_ate, function(df) df$model_ate); 
  sg_D_y_te_ate_un <- lapply(sg_D_y_te_ate, function(df) df$model_ate); 
  
  sg_D_y_tr_att_un <- lapply(sg_D_y_tr_att, function(df) df$model_att); 
  sg_D_y_te_att_un <- lapply(sg_D_y_te_att, function(df) df$model_att); 
  
  sg_D_y_tr_cru_un <- lapply(sg_D_y_tr_cru, function(df) df$model_cru);
  sg_D_y_te_cru_un <- lapply(sg_D_y_te_cru, function(df) df$model_cru); 
  
  #combine glm model from each subgroup into 1 df
  D_tr_ate <-bind_rows(sg_D_y_tr_ate_un, .id = "column_label");   
  D_te_ate <-bind_rows(sg_D_y_te_ate_un, .id = "column_label");   
  
  D_tr_att <-bind_rows(sg_D_y_tr_att_un, .id = "column_label");   
  D_te_att <-bind_rows(sg_D_y_te_att_un, .id = "column_label");   
  
  D_tr_cru <-bind_rows(sg_D_y_tr_cru_un, .id = "column_label");   
  D_te_cru <-bind_rows(sg_D_y_te_cru_un, .id = "column_label");   
  
  all_list <- list(Dlt_Y=Dlt_Y,
                   D_tr_ate = D_tr_ate,  
                   D_te_ate = D_te_ate,  
                   
                   D_tr_att = D_tr_att, 
                   D_te_att = D_te_att,  
                   
                   D_tr_cru = D_tr_cru,  
                   D_te_cru = D_te_cru  )
  
  return(all_list)
}






#' MAKE_GROUP_COLUMN
#' 
#' Function that make intearction term of W:subgroups(dummy varialbe)  
#' @param dat the dataset with subgroup ID and defintion for each observation
#' @param depth the label describe type (depth) of iCF dataset
#' @param dummytype the column created could be W:group, group, or group&W:group
#' @param level if level="g", then make dummy variable for each specific subgroup g; if level="G" then rename "SubgroupID" as "G_i"
#' 
#' @return the dataset with interaction terms of W:subgroups
#' @export
#' 
MAKE_GROUP_COLUMN <- function(dat, depth,dummytype, level){
  
  if (level=="g"){
    
    for(i in 2:(length(levels(as.factor(dat$SubgroupID))) ) ){
      
      
      if (dummytype=="W:group") {
        dat <- dat %>% dplyr:: mutate(!!paste0("W_", depth, "_g_", i)  := ifelse( SubgroupID == i, 1*W, 0*W )) 
      } else if (dummytype=="group"){      
        dat <- dat %>%     dplyr:: mutate(!!paste0(depth, "_g_", i, ": ", levels(as.factor(dat$Definition))[i])  := structure( ifelse( SubgroupID == i, 1, 0 ), label= levels(as.factor(dat$Definition))[i] )  ) 
      } else if (dummytype=="group&W:group"){
        dat <- dat %>%     dplyr:: mutate(!!paste0(depth, "_g_", i)  := ifelse( SubgroupID == i, 1, 0 ))  %>% dplyr:: mutate(!!paste0("W_", depth, "_group_", i)  := ifelse( SubgroupID == i, 1*W, 0*W )) 
        
      }
    }
    
  } 
  
  #rename SubgroupID, Defintion to prepare for merge
  dat2 <- dat %>% dplyr::mutate(Definition = as.factor(Definition),
                                SubgroupID = as.factor(SubgroupID)) %>%
    dplyr::rename(!!paste0("Definition","_", depth)  := Definition,
                  !!paste0("SubgroupID","_", depth)  := SubgroupID) %>%
    dplyr::arrange((ID))
  
  return(dat2)
}


#' MAKE_GROUP_COLUMN
#' 
#' Function that make intearction term of W:subgroups(dummy varialbe)  
#' @param dat the dataset with subgroup ID and defintion for each observation
#' @param depth the label describe type (depth) of iCF dataset
#' @param dummytype the column created could be W:group, group, or group&W:group
#' @param level if level="g", then make dummy variable for each specific subgroup g; if level="G" then rename "SubgroupID" as "G_i"
#' 
#' @return the dataset with interaction terms of W:subgroups
#' @export
#' 
MAKE_GROUP_COLUMN <- function(dat, depth,dummytype, level){
  
  if (level=="g"){
    
    for(i in 2:(length(levels(as.factor(dat$SubgroupID))) ) ){
      
      
      if (dummytype=="W:group") {
        dat <- dat %>% dplyr:: mutate(!!paste0("W_", depth, "_g_", i)  := ifelse( SubgroupID == i, 1*W, 0*W )) 
      } else if (dummytype=="group"){      
        dat <- dat %>%     dplyr:: mutate(!!paste0(depth, "_g_", i, ": ", levels(as.factor(dat$Definition))[i])  := structure( ifelse( SubgroupID == i, 1, 0 ), label= levels(as.factor(dat$Definition))[i] )  ) 
      } else if (dummytype=="group&W:group"){
        dat <- dat %>%     dplyr:: mutate(!!paste0(depth, "_g_", i)  := ifelse( SubgroupID == i, 1, 0 ))  %>% dplyr:: mutate(!!paste0("W_", depth, "_group_", i)  := ifelse( SubgroupID == i, 1*W, 0*W )) 
        
      }
    }
    
  } 
  
  #rename SubgroupID, Defintion to prepare for merge
  dat2 <- dat %>% dplyr::mutate(Definition = as.factor(Definition),
                                SubgroupID = as.factor(SubgroupID)) %>%
    dplyr::rename(!!paste0("Definition","_", depth)  := Definition,
                  !!paste0("SubgroupID","_", depth)  := SubgroupID) %>%
    dplyr::arrange((ID))
  
  return(dat2)
}


STATTEST <- function(dat ){
  dat %>% mutate(SubgroupID=as.factor(SubgroupID))
  
  X_cell = model.matrix(~ factor(SubgroupID)-1, data = dat)
  fit = lm(Y ~ factor(SubgroupID)-1, data = dat) 
  summary(fit)
  #Levene's Test for Homogeneity of Variance 
  HOV <- car::leveneTest(y=dat$Y, group=dat$SubgroupID, SubgroupID=median)
  
  #One-way analysis of means (not assuming equal variances)
  OWAM <- oneway.test(Y ~ SubgroupID, data= dat, var.equal=FALSE)
  
  #Posthoc multiple comparisons of means: Scheffe Test 
  Pairwise <- DescTools::ScheffeTest(aov(Y ~ factor(SubgroupID), data = dat))
  Pairwise_df <- as.data.frame(Pairwise$`factor(SubgroupID)` )
  
  #Bonferroni corrected P-value cutoff
  P_cutoff = 0.05/length(unique(dat$SubgroupID))
  
  # % of P-values below cutoff (i.e. pass cutoff)
  Pct_pass_cutoff <- nrow( Pairwise_df %>% dplyr::filter(pval <= P_cutoff ) )/nrow(Pairwise_df)
  
  return(list(HOV= HOV, OWAM  = OWAM, Pairwise=Pairwise, P_cutoff=P_cutoff, Pct_pass_cutoff= Pct_pass_cutoff))
}



#' Function that combine Delt_Y and its confidence interval, 2 sample test (W=1 vs W=0) in each subgroups
#' @param DltY_iCF list containing ATT, ATE of D4,3,2 subgroups in training and testing data
#' 
#' @return the dataset that combines Delt_Y and its confidence interval, 2 sample test (W=1 vs W=0) in each subgroups
#' @export
#' 
DLTYCOMBO <- function(DltY_iCF){
  
  Dlt_Y_iCF_df    <- DltY_iCF$Dlt_Y; 
  Dlt_Y_List <- list (Dlt_Y_iCF_df    = Dlt_Y_iCF_df)
  Dlt_Y_List1 <- lapply(Dlt_Y_List, function(df) nrow(df))
  nrow_max <- max(unlist(Dlt_Y_List1))
  #repeat the last row of Dlt_Y_IT_df the as many times as the row # difference with the df with largest row #  
  Dlt_Y_List2 <- lapply(Dlt_Y_List, function(df) rbind(df, transform(df[rep(nrow(df), abs(nrow_max- nrow(df) )  ),] ) ) )
  DltY_iCF_df2    <- Dlt_Y_List2$Dlt_Y_iCF_df  
  
  return(DltY_iCF_df2)
}



get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

REMOVE_LEGEND <- function(input_gg) {
  output_gg    <- 	input_gg     + theme(legend.position="none")
  return(output_gg)
}


GG_MERGE_2h <- function(fig1,fig2,Name){
  tiff(paste(Name, treeNo, iterationNo, Ntrain, intTRUE, ".tiff", sep="_"), units="in",width=15, height=10, res=150)
  g4 <- cowplot::plot_grid( fig1, fig2,
                            ncol  = 2, labels = "AUTO", label_size = 20
  )
  dev.off()
  return(g4)
}



#---------------------------------------------------------------
#  SPLIT FREQUENCY HEAT MAP
#---------------------------------------------------------------



GG_SPLITFREQ_HEAT <- function(iCF_D_SF_df, tree_depth, type_i, iflegend, truth_description, vars, iterationNo){
  
  iCF_D_SF_df_NOK <- lapply(iCF_D_SF_df, function(df) df[ , !(names(df) %in% c("k"))])
  #distinguish between iteration=1, ALL, and majority iterations
  Stability <-  format(round(length(iCF_D_SF_df)/iterationNo, digits=2), nsmall = 2)    
  #-----------------------------------------
  # to distinguish raw, iCF, and voted iCF
  #-----------------------------------------
  if (type_i=="i1"){
    #only 1 dataframe, no iteration!!!
    D = 	iCF_D_SF_df_NOK
    iteration_char = ""
    CF = "CF"
    #    subTITLE = paste0(#"No. of forest: ", length(iCF_D_SF_df) 
    #                      "Split Frequency"
    #                      ) 
  } else if(type_i=="all") {
    #get mean of MATRIX, c(2,3) means transforms the data to a 3-dimensional array and then takes a column mean out of it
    D = plyr::aaply(plyr::laply(iCF_D_SF_df_NOK, as.matrix), c(2, 3), mean)	
    iteration_char = ""
    CF = "iCF"
    #    subTITLE = paste0("No. of forests: ", length(iCF_D_SF_df) ) 
  } else if(type_i=="majority" & length(iCF_D_SF_df) == 1){
    D = 	iCF_D_SF_df_NOK
    iteration_char = "from voted best tree"
    CF= "iCF"
    #    subTITLE = paste0("No. of forest: ", length(iCF_D_SF_df), "; stability score = ",  Stability)  
  } else {
    D = plyr::aaply(plyr::laply(iCF_D_SF_df_NOK, as.matrix), c(2, 3), mean)	
    iteration_char = "from voted best tree"
    CF="iCF"
    #    subTITLE = paste0("No. of forests: ", length(iCF_D_SF_df), "; stability score = ",  Stability )  
  }
  
  if(tree_depth=="D2"){
    breaks_para <- c(-2,-1)
    title_para = ""
    subtitle_para = paste(tree_depth, CF, iteration_char)
    size_para=5
    label_para= c("2","1")
    
  } else if (tree_depth=="D3") {
    breaks_para <- c(-3,-2,-1)
    title_para = ""
    subtitle_para = paste(tree_depth, CF, iteration_char)
    size_para=5
    label_para= c("3","2","1")
    
    
  } else if (tree_depth=="D4") {
    breaks_para <- c(-4,-3, -2, -1)
    title_para = ""
    subtitle_para = paste(tree_depth, CF, iteration_char)
    size_para=5
    label_para= c("4", "3","2","1")
    
  } else {
    breaks_para <- c(-5, -4,-3, -2, -1)
    title_para = paste(truth_description)
    subtitle_para = "raw CF"
    size_para=4
    label_para= c("5", "4", "3","2","1")
    
  }
  
  if(iflegend=="nolegend"){
    legend_para = "none"} else {
      legend_para = "right"}
  
  
  d <- data.frame(D)
  p <- length(vars)
  max.depth_para= as.numeric(stringr::str_sub(tree_depth,-1,-1))
  
  #wide to long for MATRIX!
  dm <- reshape2::melt(d, id.vars = names(dimnames(d)), value.name = "value")
  
  dm$depth = rep(1:max.depth_para, p)
  
  #the line below is wrong code because "sort" function changed the variable order, the split frequency doesn't match variables anymore!!!
  #wrong code: dm <- data.frame(variable = sort(rep(names(d), nrow(d))), value = as.vector(as.matrix(d)), depth = rep(1:max.depth_para, p))
  
  for(i in 1:max.depth_para){
    tot.depth <- sum(dm[dm$depth == i,]$value)
    dm[dm$depth == i,]$value <- dm[dm$depth == i,]$value / tot.depth
  }
  
  
  g2 <- ggplot(dm, aes(x=variable, 
                       y=-depth, 
                       fill = value)) +
    geom_tile() +
    xlab("Variable") +
    ylab("Depth") +
    viridis::scale_fill_viridis(limits=c(0,1),discrete=FALSE) + 
    ggtitle( title_para,
             #paste(#tree_depth, CF, iteration_char  
             #truth_description
             #), 
             subtitle =   subtitle_para #paste(tree_depth, CF, iteration_char)
    ) + 
    theme(
      panel.background = element_blank(),
      plot.title    = element_text(color = "black", size=18, face="bold", # hjust = 0.5, 
                                   vjust=5, 
                                   margin =ggplot2:: margin(0.5, 0, 0, 0, "cm")),
      plot.subtitle = element_text(color = "black", size=16, face="bold" ),
      axis.title.x  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
      axis.title.y  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
      plot.caption  = element_text(color = "black", size=12, face="plain"),
      axis.text.x   = element_text(color = "black", size=14, face="plain", angle = 0),
      axis.text.y   = element_text(color = "black", size=14, face="plain", angle = 0)) +
    scale_y_continuous(breaks=breaks_para, expand = c(0, 0), labels = label_para) +  
    scale_x_discrete(expand = c(0, 0)) +
    theme(legend.position = legend_para) +
    geom_text(aes(label = round(value, 2)), size = size_para) 
  
  return(g2)
  
}


#---------------------------------------------------------------
#Generating figure summary of HTE P-value for D4,D3,D2 iCF
#---------------------------------------------------------------

#' Function that produce P-value for miCF
#' @param df_D4 depth 4 CF 
#' @param df_D3 depth 3 CF
#' @param df_D2 depth 2 CF
#' @param type_i from all raw CF ("all") or those CF with best trees with the same structure as those majority voted best trees ("majority")
#' @param truth_describe describe the truth
#' @param iflegend if add legend
#' 
#' @return the dataframe of residual & MSE
#' 
#' @export
#' 
#Example: all_P_HTE   <- GG_HTE_P(HTE_P_iCF_D4,   HTE_P_iCF_D3,   HTE_P_iCF_D2,   "all",    truth_description, "nolegend")



#lab over

GG_HTE_P <- function(df_D4, df_D3, df_D2, type_i, truth_describe, iflegend){
  #stack P-value of different depth of iCF
  HTE_P_iCF_D432 <-rbind(df_D4, df_D3, df_D2)
  #make a factor to make sure y-aixs shown as D4 then D3 then D2
  HTE_P_iCF_D432$tree_depth_f <- factor(HTE_P_iCF_D432$tree_depth, levels =c("D4","D3","D2"), labels = c("D4 iCF","D3 iCF","D2 iCF"))
  #obtain depth of iCF with the maxium occurence (i.e. Stability) 
  df_count<- HTE_P_iCF_D432 %>% dplyr::count(tree_depth_f) %>% slice(which.max(n)) #maybe reload dplyr as it is affected by plyr
  HTE_P_iCF_D432$rawP <- HTE_P_cf.raw
  
  if (type_i == "majority") {
    title_para=paste0(#"iCF from voted best tree for ",
      truth_describe)
    subTITLE = paste0("No. of forests for D4, D3, and D2 iCF: ", nrow(df_D4),", ", nrow(df_D3),", and ", nrow(df_D2), ", respectively" ) 
  } else  {
    title_para=paste0(#"iCF for ",
      truth_describe)
    subTITLE = paste0("No. of forests: ", nrow(df_D4) ) 
    
  }
  
  scale_para = df_count$n/iterationNo
  #if (df_count$n/iterationNo < 0.15 ) {
  #  scale_para=0.1
  #} else  {
  #  scale_para=1
  #}
  
  
  if(iflegend=="nolegend"){
    legend_para = "none"}
  else {legend_para = "right"}
  
  if(mean( round(  HTE_P_iCF_D432$rawP) )<0.001){
    Pvalue_para = "< 0.001"
    caption_para = "The raw CF P-value "
  } else {
    Pvalue_para = mean( round(  HTE_P_iCF_D432$rawP ,3) )
    caption_para = "The raw CF P-value = "
    
  }
  
  theme_set(theme_minimal())
  P_ridge <- ggplot(HTE_P_iCF_D432[ which(HTE_P_iCF_D432$HTE_P_cf >= 0), ],  
                    aes(x =HTE_P_cf, y=forcats::fct_rev(tree_depth_f), fill = stat(x)) ) +
    geom_density_ridges_gradient(scale =scale_para, #A scaling factor to scale the height of the ridgelines relative to the spacing between them. 
                                 #A value of 1 indicates that the maximum point of any ridgeline touches the baseline right above, 
                                 #assuming even spacing between baselines.
                                 size = 0.3, 
                                 rel_min_height = 0 #Lines with heights below this cutoff will be removed. 
    ) +
    scale_fill_viridis_c(name = "P-value", option = "D", begin = 0, end = 1,limits=c(0,1) ) +
    labs(title = title_para, 
         subtitle = , #subTITLE ,
         x="P-value", 
         y="Density" ,
         caption = paste0(caption_para, Pvalue_para, "(dotted line)") 
    ) +
    theme(
      plot.title    = element_text(color="black", size=18, face="bold"),
      plot.subtitle = element_text(color="black", size=16, face="bold"),
      axis.title.x  = element_text(color="black", size=14, face="plain"),
      axis.title.y  = element_text(color="black", size=14, face="plain"),
      plot.caption  = element_text(color="black", size=12, face="plain"),
      axis.text.x   = element_text(color="black", size=14, face="plain", angle=0),
      axis.text.y   = element_text(color="black", size=14, face="plain", angle=0))+
    #don't activate: geom_vline(aes(xintercept =   HTE_P_iCF_D432$rawP ),col='gray45', size=0.6 ,linetype="dashed") + #doesn't work when merging ggplot as the lastest vline will be applied to all other ggplots
    theme(legend.position = legend_para)
  return(P_ridge)
}


PlotTrimHDcov <- function(Accuracy_X){
  VS_label <- names(Accuracy_X)
  plot(unlist(Accuracy_X), main="The effect of trimming HD covariates on \n accuracy of predicting treatment assignment",
       xlab=expression(bold("Cutoff percentile of variable importance value")), 
       ylab=expression(bold("Accuracy")),
       xaxt='n',
       lwd=3, 
       type='l'
  )
  points(unlist(Accuracy_X), bg='tomato2', pch=21, cex=2, lwd=1) #
  axis(1, at=1:19, labels=VS_label)
  
}



#' Function that shows subgroup Delta Y by subgrooup decisions obtained from each depth of forest (oneCF, iCF, or iCFv) 
#' @param model_df the df of coefficient (i.e. Delta Y) and its CI  
#' @methodL forest type: oneCF, iCF, or iCFv
#' @legend_para whether show legend of gradient fill, "nolgend" or "right"
#' @splitter splitter type
#' 
#' @return The ggplot
#' 
#' @export
#' 
GG_SG_DELTA_Y <- function(model_df, tree_depth, legend_para, splitter){
  
  if (model_df$Sub_def[1]=="W <100"){model_df$Sub_def="No Heterogeneous Subgroups"}
  
  if (splitter=="allBinary"){
    #first deal with the non-last subgroup definition
    model_df <- model_df %>% dplyr::mutate(Sub_def =   gsub ("<=0 &", "=0 &",Sub_def) )%>%
      dplyr::mutate(Sub_def =   gsub (">0 &", "=1 &", Sub_def) )
    Sub_def_L <- as.list(model_df$Sub_def) 
    #str_sub doesn't work in column of dataframe or a list, thus have to use paste0 Fx to combine first and second part
    Sub_def_L_last_3 <- lapply(Sub_def_L, function(E) if(str_sub(E, -3, -1)=="<=0"){paste0( str_sub(E,1,-4), stringr::str_sub(E, -3,-1) <- "=0" ) } else{E} )
    
    Sub_def_L_last_3_2 <- lapply(Sub_def_L_last_3, function(E) if(str_sub(E, -2, -1)==">0"){paste0( str_sub(E,1,-3), stringr::str_sub(E, -2,-1) <- "=1" ) } else{E} )
    
    Sub_def_v = unlist(Sub_def_L_last_3_2, use.names=FALSE) 
    
    model_df <- model_df %>% mutate(Sub_def_v = Sub_def_v)
    
  }
  
  if(tree_depth=="D2" | tree_depth=="D3" | tree_depth=="D4"){
    title_para = ""
  }  else if (tree_depth=="truth") {
    title_para = paste(truth_description)
    
  } 
  
  g3<- ggplot(model_df, aes(Estimate, column_label, color= )) + 
    geom_point(aes(x=Estimate, color = Estimate), size=5) + 
    geom_linerange(aes(xmin = conf.low, xmax=conf.high, color = Estimate),   size=1) + 
    geom_text(aes(label=Sub_def_v), vjust = 0, nudge_y = 0.2, size=2.5) +
    scale_color_viridis(option = "D", name = expression(paste(Delta, "Y")))+
    ggtitle(title_para) + 
    labs(subtitle= ifelse(tree_depth=="truth", paste("Truth"), paste(tree_depth, "iCF", sep = " ")), 
         caption = ifelse(tree_depth=="truth", paste(""),      paste0("Stability = ", mean(model_df$stability) ))  
         
    ) +
    ylab("Subgroup") + xlab(expression(paste("Subgroup ", Delta, "Y"))) +
    theme_bw() +
    theme(
      plot.title    = element_text(color = "black", size=16.5, face="bold", # hjust = 0.5, 
                                   vjust=5, 
                                   margin =ggplot2:: margin(0.5, 0, 0, 0, "cm")),
      plot.subtitle = element_text(color="black", size=16, face="bold"),
      axis.title.x  = element_text(color="black", size=14, face="plain"),
      axis.title.y  = element_text(color="black", size=14, face="plain"),
      plot.caption  = element_text(color="black", size=14, face="plain"),
      axis.text.x   = element_text(color="black", size=14, face="plain", angle=0),
      axis.text.y   = element_text(color="black", size=14, face="plain", angle=0))  + 
    scale_y_discrete(limits = unique(rev(model_df$column_label))) + 
    theme(legend.position = legend_para #"right"
          #"nolegend"
    )
  return(g3)
}


#'PLOT_BT
#'
#' This function provides the plot of best tree
#' @param forest causal forest
#' @param type "causal" or "instrumental"
#' 
#' @return a list of tree info and best tree plot
#'
#' @export
#'

PLOT_BT<-function(forest, type){
  best_tree_info<-find_best_tree(forest, type)
  best_tree_info$best_tree
  plot_bt<-plot(grf::get_tree(forest, best_tree_info$best_tree))
  return(list(best_tree_info$best_tree, plot_bt))
}


#' ROC2comp
#' 
#' Function that compares 2 ROCs, not necessary for iCF/hdiCF
#' @param response1 true outcome 1
#' @param response2 true outcome 2
#' @param predict1 prediction 1
#' @param predict2 prediction 2  
#' @param label1 label 1
#' @param label2 lebel 2 
#' @param outcome outcome, e.g. Y or W
#' @return the image object
#' 
#' @export


ROC2comp <- function(response1, predict1, response2, predict2,label1, label2, outcome){
  dev.off()
  
  par(pty = "s") 
  roc_Y.hat_all <-pROC::roc(response = response1, predictor = predict1, plot=T, 
                            legacy.axes=TRUE, percent=TRUE, 
                            xlab="False Positive Percentage", ylab="True Postive Percentage", 
                            col="#377eb8", lwd=4, print.auc=TRUE)
  
  pROC::plot.roc(response2, predict2, percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
  legend("bottomright", 
         legend=c(label1, label2), 
         col=c("#377eb8", "#4daf4a"), 
         lwd=4)
  title(main =  paste0("ROC of predicted ", outcome) ,  line = 3, adj = 0)
  
  ROC_compare <- recordPlot()
  return( ROC_compare )
  
}
