#############################################################################################
# iCF CV for entire process
# Author: Tiansheng Wang  
# Last update date:3/1/2022
# Version: 0.1         
#############################################################################################

SUBGROUP_PIPLINE<- function(leafsize, treeNo, iterationNo, Ntrain, Train_ID, Test_ID ){
  
  s_iCF=Sys.time()
  
  if (iterationNo>1){
    iCF_D4<- iCF(leafsize$D4,  treeNo, iterationNo, Ntrain, "D4",  split_val_round_posi)
    iCF_D3<- iCF(leafsize$D3,  treeNo, iterationNo, Ntrain, "D3",  split_val_round_posi) 
    #when leaf size too small, D2 has this error "Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) : no rows to aggregate"
    #thus to make it auto run, if iCF_D2 has error, then make it equal to iCF_D3 temparily.
    if ( is.null( tryCatch( iCF(leafsize$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi), error=function(e){}) ) == TRUE ) {
      warning("The minimum leaf size is too small for Depth 2 causal forest, need to increase it!")
      iCF_D2 = iCF_D3
    } else {
      iCF_D2 = iCF(leafsize$D2,  treeNo, iterationNo, Ntrain, "D2",  split_val_round_posi)
    }
    
    #best tree dataframe	
    iCF_D4_BT <<- GET_TREE_L(iCF_D4, 1) 
    iCF_D3_BT <<- GET_TREE_L(iCF_D3, 1)
    iCF_D2_BT <<- GET_TREE_L(iCF_D2, 1) 
    
    #mean number of node from a causal tree dataframe
    mean_node_D4 <-round(mean( do.call("rbind" , lapply(iCF_D4_BT, function(df) nrow(df))) ), 0)
    mean_node_D3 <-round(mean( do.call("rbind" , lapply(iCF_D3_BT, function(df) nrow(df))) ), 0)
    mean_node_D2 <-round(mean( do.call("rbind" , lapply(iCF_D2_BT, function(df) nrow(df))) ), 0)
    
    
    #best tree list format
    iCF_D4_BT_L <- GET_TREE_L(iCF_D4, 2) 
    iCF_D3_BT_L <- GET_TREE_L(iCF_D3, 2)
    iCF_D2_BT_L <- GET_TREE_L(iCF_D2, 2)
    
    
    #for split frequency, doesn't consider split values		
    iCF_D4_SF <- GET_TREE_L(iCF_D4, 3)
    iCF_D3_SF <- GET_TREE_L(iCF_D3, 3)
    iCF_D2_SF <- GET_TREE_L(iCF_D2, 3)
    
    #HTE P value for ALL iteractions of CF
    HTE_P_iCF_D4 <- HTE_P_extract(iCF_D4_BT, "D4")
    HTE_P_iCF_D3 <- HTE_P_extract(iCF_D3_BT, "D3")
    HTE_P_iCF_D2 <- HTE_P_extract(iCF_D2_BT, "D2")
    
    HTE_P_D4_median <- 	median(HTE_P_iCF_D4$HTE_P_cf)
    HTE_P_D3_median <- 	median(HTE_P_iCF_D3$HTE_P_cf)
    HTE_P_D2_median <- 	median(HTE_P_iCF_D2$HTE_P_cf)
    
  } else if (iterationNo==1){
    
    iCF_D4 <- oneCF(leafsize$D4,  treeNo, "none", "D4", split_val_round_posi) 
    iCF_D3 <- oneCF(leafsize$D3,  treeNo, "none", "D3", split_val_round_posi) 
    iCF_D2 <- oneCF(leafsize$D2,  treeNo, "none", "D2", split_val_round_posi)  
    
    #tree b dataframe	
    iCF_D4_BT <- iCF_D4$besttreelist              
    iCF_D3_BT <- iCF_D3$besttreelist              
    iCF_D2_BT <- iCF_D2$besttreelist 
    #for split frequency		
    iCF_D4_SF <- list(iCF_D4$d)          
    iCF_D3_SF <- list(iCF_D3$d)           
    iCF_D2_SF <- list(iCF_D2$d)
    
  }
  
  
  iCF_D4_t  <-PRE_MAJORITY_TREE  (iCF_D4_BT)
  iCF_D3_t  <-PRE_MAJORITY_TREE  (iCF_D3_BT)
  iCF_D2_t  <-PRE_MAJORITY_TREE  (iCF_D2_BT)
  
  iCF_D4_t_r<-PRE_MAJORITY_TREE_R(iCF_D4_BT)
  iCF_D3_t_r<-PRE_MAJORITY_TREE_R(iCF_D3_BT)
  iCF_D2_t_r<-PRE_MAJORITY_TREE_R(iCF_D2_BT)
  
  vote_D4_tree       <-MAJORITY_VOTE(iCF_D4_BT,         iCF_D4_t,        iCF_D4_t,        iCF_D4_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  vote_D3_tree       <-MAJORITY_VOTE(iCF_D3_BT,         iCF_D3_t,        iCF_D3_t,        iCF_D3_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  vote_D2_tree       <-MAJORITY_VOTE(iCF_D2_BT,         iCF_D2_t,        iCF_D2_t,        iCF_D2_SF,      tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
  
  vote_D4_tree_R  <-MAJORITY_VOTE(iCF_D4_BT,      iCF_D4_t_r,   iCF_D4_t_r,   iCF_D4_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  vote_D3_tree_R  <-MAJORITY_VOTE(iCF_D3_BT,      iCF_D3_t_r,   iCF_D3_t_r,   iCF_D3_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  vote_D2_tree_R  <-MAJORITY_VOTE(iCF_D2_BT,      iCF_D2_t_r,   iCF_D2_t_r,   iCF_D2_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
  
  vote_D4_tree.syn     = vote_D4_tree_R$majority.syn
  vote_D3_tree.syn.ori = vote_D3_tree_R$majority.syn
  vote_D2_tree.syn.ori = vote_D2_tree_R$majority.syn
  
  HTE_P_D4_M_median <- vote_D4_tree$HTE_P_cf_M_median
  HTE_P_D3_M_median <- vote_D3_tree$HTE_P_cf_M_median
  HTE_P_D2_M_median <- vote_D2_tree$HTE_P_cf_M_median
  
  ###############################################
  # original split value without imputation
  ###############################################
  vote_D4_subgroup <- TREE2SUBGROUP(vote_D4_tree.syn )
  vote_D3_subgroup <- TREE2SUBGROUP(vote_D3_tree.syn.ori )
  vote_D2_subgroup <- TREE2SUBGROUP(vote_D2_tree.syn.ori )
  
  
  
  
  s_SGdeci=Sys.time()
  
  DATA4grplasso_iCF_tr     <- DATA4grplasso(vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup, Train_ID)
  DATA4grplasso_iCF_te     <- DATA4grplasso(vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup, Test_ID)
  
  Train_ID_SG_iCF       <- DATA4grplasso_iCF_tr$Dat_ID_SG
  Test_ID_SG_iCF        <- DATA4grplasso_iCF_te$Dat_ID_SG
  
  #the following formulas work for both binary and continuous outcome
  #actually, formula_gl = formula_g234; and iCF and oneCF can use the same formula, thus did not distinguish the name
  formula_basic<<-  as.formula(  paste0("Y ~ W +", paste0(colnames(X), collapse = " + ") )   )
  formula_gl   <<- DATA4grplasso_iCF_tr$formula_gl
  formula_g2   <<- DATA4grplasso_iCF_tr$formula_g2
  formula_g3   <<- DATA4grplasso_iCF_tr$formula_g3
  formula_g4   <<- DATA4grplasso_iCF_tr$formula_g4
  formula_g23  <<- DATA4grplasso_iCF_tr$formula_g23
  #_________________________________________________
  
  
  Deci_Final_iCF.act.tr        <- CF_GROUP_DECISION(HTE_P_cf.raw, Train_ID_SG_iCF, "iCF", "actual", vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup)
  Deci_Final_iCF.tran.tr       <- CF_GROUP_DECISION(HTE_P_cf.raw, Train_ID_SG_iCF, "iCF", "transform", vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup)
  Deci_Final_iCF.act.te        <- CF_GROUP_DECISION(HTE_P_cf.raw, Test_ID_SG_iCF, "iCF", "actual", vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup)
  Deci_Final_iCF.tran.te       <- CF_GROUP_DECISION(HTE_P_cf.raw, Test_ID_SG_iCF, "iCF", "transform", vote_D4_subgroup, vote_D3_subgroup, vote_D2_subgroup)
  
  
  
  e_iCF=Sys.time()
  
  
  
  return(list ( vote_D4_tree.syn = vote_D4_tree.syn,
                vote_D3_tree.syn = vote_D3_tree.syn.ori,
                vote_D2_tree.syn =  vote_D2_tree.syn.ori,
                vote_D4_subgroup=vote_D4_subgroup,
                vote_D3_subgroup=vote_D3_subgroup,
                vote_D2_subgroup=vote_D2_subgroup,
                time_iCF_AIC =  difftime(e_iCF, s_iCF, units="secs"), 
   
           
                stability_D4_T       = unlist(vote_D4_tree$stability),        
                stability_D4_T_r     = unlist(vote_D4_tree_R$stability),      
                stability_D3_T       = unlist(vote_D3_tree$stability),        
                stability_D3_T_r     = unlist(vote_D3_tree_R$stability),       
                stability_D2_T       = unlist(vote_D2_tree$stability),   
                stability_D2_T_r     = unlist(vote_D2_tree_R$stability),
                Deci_Final_iCF.act.tr   = Deci_Final_iCF.act.tr,
                Deci_Final_iCF.tran.tr  = Deci_Final_iCF.tran.tr ,
                Deci_Final_iCF.act.te  = Deci_Final_iCF.act.te,
                Deci_Final_iCF.tran.te  = Deci_Final_iCF.tran.te ,
                grp_AIC_pick_D_iCF.act.tr     = Deci_Final_iCF.act.tr$grp_AIC_pick_D,
                grp_AIC_pick_D_iCF.tran.tr    = Deci_Final_iCF.tran.tr$grp_AIC_pick_D,
                grp_AIC_pick_D_iCF.act.te     = Deci_Final_iCF.act.te$grp_AIC_pick_D,
                grp_AIC_pick_D_iCF.tran.te    = Deci_Final_iCF.tran.te$grp_AIC_pick_D,
                Deci_AIC_iCF.act.tr     = Deci_Final_iCF.act.tr$Deci_Final_CF_AIC,
                Deci_AIC_iCF.tran.tr    = Deci_Final_iCF.tran.tr$Deci_Final_CF_AIC,
                Deci_AIC_iCF.act.te     = Deci_Final_iCF.act.te$Deci_Final_CF_AIC,
                Deci_AIC_iCF.tran.te    = Deci_Final_iCF.tran.te$Deci_Final_CF_AIC,
                #Train_ID_SG_iCF      = Train_ID_SG_iCF,
                Test_ID_SG_iCF       = Test_ID_SG_iCF 
                
  )
  )  
}




TREE_PERFORMANCE <- function(Deci_SG, iterationNo){
  #performance 
  DltY_iCF_te      = DltY_DATA(Deci_SG,       Test_ID, "predict" )
  MSE_iCF_te       <- MSE_DltY(DltY_true_te, DltY_iCF_te)
  MSE_ate_iCF_te   <- as.numeric(MSE_iCF_te$MSE_ate)
  MSE_att_iCF_te   <- as.numeric(MSE_iCF_te$MSE_att)
  
  #DISCOVERYRATE function take care of the "NA" decision
  if(iterationNo >1){
    disco_SG_iCF       = DISCOVERYRATE(Deci_SG ,                 truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")$value
    disco_INT_iCF      = DISCOVERYRATE(SG2INT(Deci_SG),          truth_description, tree_true_subgroup, truth_INT, "iCF", "interaction")$value
    #disco_SG_iCF.act.aa    = DISCOVERYRATE(Deci_SG,                  truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")$sg_Nacc_Nall
    #disco_SG_iCF.act.aat   = DISCOVERYRATE(Deci_SG,                  truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")$sg_Nacc_NallT
    disco_CATEmax_iCF      = DISCOVERYRATE(CATE_MAX_SG(DltY_iCF_te), truth_description, tree_true_subgroup, CATE_max_true_te, "iCF", "CATEmax")$value
    
  } else if (iterationNo==1) {
    disco_SG_iCF       = DISCOVERYRATE(Deci_SG ,                  truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$value
    disco_INT_iCF      = DISCOVERYRATE(SG2INT(Deci_SG),           truth_description, tree_true_subgroup, truth_INT, "oneCFb", "interaction")$value
    #disco_SG_iCF.act.aa    = DISCOVERYRATE(Deci_SG,                   truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$sg_Nacc_Nall
    #disco_SG_iCF.act.aat   = DISCOVERYRATE(Deci_SG,                   truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$sg_Nacc_NallT
    disco_CATEmax_iCF      = DISCOVERYRATE(CATE_MAX_SG(DltY_iCF_te),  truth_description, tree_true_subgroup, CATE_max_true_te, "oneCFb", "CATEmax")$value
    
  }  
  return(list(disco_SG_iCF   = disco_SG_iCF ,
              disco_INT_iCF  = disco_INT_iCF, 
              #disco_SG_AIC_iCF.act.aa = disco_SG_AIC_iCF.act.aa,
              #disco_SG_AIC_iCF.act.aat  = disco_SG_AIC_iCF.act.aat,
              #DltY_iCF_te       = DltY_iCF_te,
              disco_CATEmax_iCF = disco_CATEmax_iCF,
              MSE_ate_iCF_te    = MSE_ate_iCF_te,
              MSE_att_iCF_te    = MSE_att_iCF_te)
  )  
  
}

                                                        
#' Function that run the iCF algorithm with CV on the whole data
#' @param K fold #
#' @param treeNo Tree # of iCF
#' @SigLevel significance level, because iCF will only be intiated when HTE P-value <0.1, thus we set significance level=1 to run iCF anyway 
#'  
#' @return subgroup decision
#' 
#' @export
  
                                                        
                                                        
iCFCV <- function(K, treeNo, SigLevel){
  s_iCF_CV=Sys.time()
  if (round(HTE_P_cf.raw,1) <= SigLevel){
    model.g2.act <- list()
    model.g3.act <- list()
    model.g4.act <- list()
    model.m.act  <- list()
    model.g2.tran <- list()
    model.g3.tran <- list()
    model.g4.tran <- list()
    model.m.tran  <- list()
    
    cf_raw_key <- list()
    mse_per_fold.g2.act     <- list()
    mse_per_fold.g3.act     <- list()
    mse_per_fold.g4.act     <- list()
    mse_per_fold.m.act      <- list()
    mse_per_fold.g2.tran    <- list()
    mse_per_fold.g3.tran    <- list()
    mse_per_fold.g4.tran    <- list()
    mse_per_fold.m.tran     <- list()
    
    accuracy_per_fold.g2.act <- list()
    accuracy_per_fold.g3.act <- list()
    accuracy_per_fold.g4.act <- list()
    accuracy_per_fold.m.act <- list()

    Train_ID_SG_iCF <- list()
    Test_ID_SG_iCF <- list()
    
    vote_D4_tree.syn <- list()
    vote_D3_tree.syn <- list()
    vote_D2_tree.syn <- list()
    vote_D4_subgroup.L <- list()
    vote_D3_subgroup.L <- list()
    vote_D2_subgroup.L <- list()
    stability_D4_T_r <- list()
    stability_D3_T_r <- list()
    stability_D2_T_r <- list()
    
    Deci_Final_iCF.act.tr <- list()
    Deci_Final_iCF.tran.tr <- list()
    Deci_Final_iCF.act.te <- list()
    Deci_Final_iCF.tran.te <- list()
    
    test_data.act <- list()
    test_data.tran <- list()
    
    Deci_AIC_test.act  <- list()
    Deci_AIC_test.tran <- list()
    
    Test_ID_SG_iCF <- list()
    
    HTE_P_cf.raw.L <- list()
    selected_cf.idx.L <- list()
    
    accuracy_per_fold <- list()

    set.seed(20160413)
    tt_indicies <- caret::createFolds(y=Train[,1], k= K)
    leafsize <<- LEAFSIZE_tune(round(Ntrain/K*(K-1),0), treeNo, iterationNo)
    #leafsize = list(D4=140, D3=70, D2=35)
    for(f in 1:length(tt_indicies)){
      
      Train_cf <- Train[-tt_indicies[[f]],]
      Test_cf  <- Train[tt_indicies[[f]],]
      ID_cf <-1:nrow(Train)
      Train_ID_cf <- cbind(Train_cf, as.vector(ID_cf[-tt_indicies[[f]]])) %>% dplyr::rename (ID=`as.vector(ID_cf[-tt_indicies[[f]]])`) 
      Test_ID_cf  <- cbind(Test_cf,  as.vector(       tt_indicies[[f]]) ) %>% dplyr::rename (ID=`as.vector(tt_indicies[[f]])`)  
      dplyr::all_equal(Train_ID_cf, Test_ID_cf)
      #============================== raw full CF for CV training data ==============================
      cf_raw_key[[f]] <- CF_RAW_key(Train_cf)   
      X      <- cf_raw_key[[f]]$X
      Y      <- cf_raw_key[[f]]$Y
      W      <- cf_raw_key[[f]]$W
      Y.hat  <- cf_raw_key[[f]]$Y.hat
      W.hat  <- cf_raw_key[[f]]$W.hat
      varimp_cf  <- cf_raw_key[[f]]$varimp_cf
      selected_cf.idx.L[[f]]  <- cf_raw_key[[f]]$selected_cf.idx
      selected_cf.idx         <- selected_cf.idx.L[[f]] 
      HTE_P_cf.raw.L[[f]] <- cf_raw_key[[f]]$HTE_P_cf.raw
      HTE_P_cf.raw <- HTE_P_cf.raw.L[[f]] 
      #==============================     iCF for CV training data   ==============================
      #SG_D[[f]]<-SUBGROUP_PIPLINE(leafsize, treeNo, iterationNo, nrow(Train_cf), Train_ID_cf, Test_ID_cf)
      
      ##################################      SUBGROUP_PIPLINE starts      ################################################################
      s_iCF=Sys.time()
      
      if (iterationNo>1){
        iCF_D4<- iCF(leafsize$D4,  treeNo, iterationNo,  nrow(Train_cf), "D4",  split_val_round_posi)
        iCF_D3<- iCF(leafsize$D3,  treeNo, iterationNo,  nrow(Train_cf), "D3",  split_val_round_posi) 
        #when leaf size too small, D2 has this error "Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) : no rows to aggregate"
        #thus to make it auto run, if iCF_D2 has error, then make it equal to iCF_D3 temparily.
        if ( is.null( tryCatch( iCF(leafsize$D2,  treeNo, iterationNo,  nrow(Train_cf), "D2",  split_val_round_posi), error=function(e){}) ) == TRUE ) {
          warning("The minimum leaf size is too small for Depth 2 causal forest, need to increase it!")
          iCF_D2 = iCF_D3
        } else {
          iCF_D2 = iCF(leafsize$D2,  treeNo, iterationNo,  nrow(Train_cf), "D2",  split_val_round_posi)
        }
        #best tree dataframe	
        iCF_D4_BT <<- GET_TREE_L(iCF_D4, 1) 
        iCF_D3_BT <<- GET_TREE_L(iCF_D3, 1)
        iCF_D2_BT <<- GET_TREE_L(iCF_D2, 1) 
        #mean number of node from a causal tree dataframe
        mean_node_D4 <-round(mean( do.call("rbind" , lapply(iCF_D4_BT, function(df) nrow(df))) ), 0)
        mean_node_D3 <-round(mean( do.call("rbind" , lapply(iCF_D3_BT, function(df) nrow(df))) ), 0)
        mean_node_D2 <-round(mean( do.call("rbind" , lapply(iCF_D2_BT, function(df) nrow(df))) ), 0)
        #best tree list format
        iCF_D4_BT_L <- GET_TREE_L(iCF_D4, 2) 
        iCF_D3_BT_L <- GET_TREE_L(iCF_D3, 2)
        iCF_D2_BT_L <- GET_TREE_L(iCF_D2, 2)
        #for split frequency, doesn't consider split values		
        iCF_D4_SF <- GET_TREE_L(iCF_D4, 3)
        iCF_D3_SF <- GET_TREE_L(iCF_D3, 3)
        iCF_D2_SF <- GET_TREE_L(iCF_D2, 3)
        #HTE P value for ALL iteractions of CF
        HTE_P_iCF_D4 <- HTE_P_extract(iCF_D4_BT, "D4")
        HTE_P_iCF_D3 <- HTE_P_extract(iCF_D3_BT, "D3")
        HTE_P_iCF_D2 <- HTE_P_extract(iCF_D2_BT, "D2")
        HTE_P_D4_median <- 	median(HTE_P_iCF_D4$HTE_P_cf)
        HTE_P_D3_median <- 	median(HTE_P_iCF_D3$HTE_P_cf)
        HTE_P_D2_median <- 	median(HTE_P_iCF_D2$HTE_P_cf)
      } else if (iterationNo==1){
        iCF_D4 <- oneCF(leafsize$D4,  treeNo, "none", "D4", split_val_round_posi) 
        iCF_D3 <- oneCF(leafsize$D3,  treeNo, "none", "D3", split_val_round_posi) 
        iCF_D2 <- oneCF(leafsize$D2,  treeNo, "none", "D2", split_val_round_posi)  
        #tree b dataframe	
        iCF_D4_BT <- iCF_D4$besttreelist              
        iCF_D3_BT <- iCF_D3$besttreelist              
        iCF_D2_BT <- iCF_D2$besttreelist 
        #for split frequency		
        iCF_D4_SF <- list(iCF_D4$d)          
        iCF_D3_SF <- list(iCF_D3$d)           
        iCF_D2_SF <- list(iCF_D2$d)
      }
      
      iCF_D4_t  <-PRE_MAJORITY_TREE  (iCF_D4_BT)
      iCF_D3_t  <-PRE_MAJORITY_TREE  (iCF_D3_BT)
      iCF_D2_t  <-PRE_MAJORITY_TREE  (iCF_D2_BT)
      iCF_D4_t_r<-PRE_MAJORITY_TREE_R(iCF_D4_BT)
      iCF_D3_t_r<-PRE_MAJORITY_TREE_R(iCF_D3_BT)
      iCF_D2_t_r<-PRE_MAJORITY_TREE_R(iCF_D2_BT)
      vote_D4_tree    <-MAJORITY_VOTE(iCF_D4_BT,      iCF_D4_t,     iCF_D4_t,     iCF_D4_SF,   tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
      vote_D3_tree    <-MAJORITY_VOTE(iCF_D3_BT,      iCF_D3_t,     iCF_D3_t,     iCF_D3_SF,   tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
      vote_D2_tree    <-MAJORITY_VOTE(iCF_D2_BT,      iCF_D2_t,     iCF_D2_t,     iCF_D2_SF,   tree_true,          tree_true_N1,   tree_true_N123, split_val_round_posi)
      vote_D4_tree_R  <-MAJORITY_VOTE(iCF_D4_BT,      iCF_D4_t_r,   iCF_D4_t_r,   iCF_D4_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
      vote_D3_tree_R  <-MAJORITY_VOTE(iCF_D3_BT,      iCF_D3_t_r,   iCF_D3_t_r,   iCF_D3_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)
      vote_D2_tree_R  <-MAJORITY_VOTE(iCF_D2_BT,      iCF_D2_t_r,   iCF_D2_t_r,   iCF_D2_SF,   tree_true_r,        tree_true_N1_r, tree_true_N123_r, split_val_round_posi)

      vote_D4_tree.syn[[f]] <- vote_D4_tree_R$majority.syn
      vote_D3_tree.syn[[f]] <- vote_D3_tree_R$majority.syn
      vote_D2_tree.syn[[f]] <- vote_D2_tree_R$majority.syn
      vote_D4_subgroup.L[[f]] <-  TREE2SUBGROUP(vote_D4_tree.syn[[f]] )
      vote_D3_subgroup.L[[f]] <-  TREE2SUBGROUP(vote_D3_tree.syn[[f]] )
      vote_D2_subgroup.L[[f]] <-  TREE2SUBGROUP(vote_D2_tree.syn[[f]] )
      HTE_P_D4_M_median <- vote_D4_tree$HTE_P_cf_M_median
      HTE_P_D3_M_median <- vote_D3_tree$HTE_P_cf_M_median
      HTE_P_D2_M_median <- vote_D2_tree$HTE_P_cf_M_median
      ###############################################
      # original split value without imputation
      ###############################################
      s_SGdeci=Sys.time()
      DATA4grplasso_iCF_tr  <- DATA4grplasso(vote_D4_subgroup.L[[f]], vote_D3_subgroup.L[[f]], vote_D2_subgroup.L[[f]], Train_ID_cf)
      DATA4grplasso_iCF_te  <- DATA4grplasso(vote_D4_subgroup.L[[f]], vote_D3_subgroup.L[[f]], vote_D2_subgroup.L[[f]], Test_ID_cf)
      Train_ID_SG_iCF[[f]]  <- DATA4grplasso_iCF_tr$Dat_ID_SG
      Test_ID_SG_iCF [[f]]  <- DATA4grplasso_iCF_te$Dat_ID_SG
       #the following formulas work for both binary and continuous outcome
      #actually, formula_gl = formula_g234; and iCF and oneCF can use the same formula, thus did not distinguish the name
      formula_basic<<-  as.formula(  paste0("Y ~ W +", paste0(colnames(X), collapse = " + ") )   )
      formula_gl   <<- DATA4grplasso_iCF_tr$formula_gl
      formula_g2   <<- DATA4grplasso_iCF_tr$formula_g2
      formula_g3   <<- DATA4grplasso_iCF_tr$formula_g3
      formula_g4   <<- DATA4grplasso_iCF_tr$formula_g4
      formula_g23  <<- DATA4grplasso_iCF_tr$formula_g23
      #_________________________________________________
      Deci_Final_iCF.act.tr[[f]]  <- CF_GROUP_DECISION(1, HTE_P_cf.raw, Train_ID_SG_iCF[[f]], "iCF", "actual",    vote_D4_subgroup.L[[f]], vote_D3_subgroup.L[[f]],  vote_D2_subgroup.L[[f]])
      Deci_Final_iCF.tran.tr[[f]] <- CF_GROUP_DECISION(1, HTE_P_cf.raw, Train_ID_SG_iCF[[f]], "iCF", "transform", vote_D4_subgroup.L[[f]], vote_D3_subgroup.L[[f]],  vote_D2_subgroup.L[[f]])
      Deci_Final_iCF.act.te[[f]]  <- CF_GROUP_DECISION(1, HTE_P_cf.raw, Test_ID_SG_iCF[[f]],  "iCF", "actual",    vote_D4_subgroup.L[[f]], vote_D3_subgroup.L[[f]],  vote_D2_subgroup.L[[f]])
      Deci_Final_iCF.tran.te[[f]] <- CF_GROUP_DECISION(1, HTE_P_cf.raw, Test_ID_SG_iCF[[f]],  "iCF", "transform", vote_D4_subgroup.L[[f]], vote_D3_subgroup.L[[f]],  vote_D2_subgroup.L[[f]])
      print("CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
      e_iCF=Sys.time()
      time_iCF_AIC =  difftime(e_iCF, s_iCF, units="secs") 
      #########################################     SUBGROUP_PIPLINE ends     ########################################################
      stability_D4_T_r[[f]] <- unlist(vote_D4_tree_R$stability)
      print("stability_D4_T_r")
      stability_D3_T_r[[f]] <- unlist(vote_D3_tree_R$stability)
      print("stability_D3_T_r")
      stability_D2_T_r[[f]] <- unlist(vote_D2_tree_R$stability)
      print("stability_D2_T_r")
      model.g2.act[[f]] =  Deci_Final_iCF.act.tr[[f]]$model.g2
      print("model.g2.act")
      model.g3.act[[f]] =  Deci_Final_iCF.act.tr[[f]]$model.g3
      print("model.g3.act")
      model.g4.act[[f]] =  Deci_Final_iCF.act.tr[[f]]$model.g4
      print("model.g4.act")
      model.m.act[[f]]  = Deci_Final_iCF.act.tr[[f]]$model.basic
      print("model.m.act")
      model.g2.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g2
      print("model.g2.tran")
      model.g3.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g3
      print("model.g3.tran")
      model.g4.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g4
      print("model.g4.tran")
      model.m.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.basic
      print("model.m.tran")
      test_data.act[[f]] = SGMODEL_DATA(Test_ID_SG_iCF[[f]], "actual")$dat_ID_SG_df
      print(" test_data.act")
      test_data.tran[[f]] =SGMODEL_DATA(Test_ID_SG_iCF[[f]], "transform")$dat_ID_SG_df
      print(" test_data.tran")
      #====================using AIC value from models on testing data====================
      Deci_AIC_test.act[[f]]  <- Deci_Final_iCF.act.te[[f]]$Deci_Final_CF_AIC
      print("Deci_AIC_test.act")
      Deci_AIC_test.tran[[f]] <- Deci_Final_iCF.tran.te[[f]]$Deci_Final_CF_AIC
      print("Deci_AIC_test.tran")
      Ncol_g234= ncol(test_data.act[[f]] %>% dplyr::select((contains( c("G4", "G3", "G2") ) ) ))      
      print("Ncol_g234")

      if (length(unique(Y)) >=8){
        test_data.act[[f]]$predict.g2.act <- predict(model.g2.act[[f]], newdata=test_data.act[[f]][,1:(2 + ncol(X) + Ncol_g234)])
        print("test_data.act[[f]]$predict.g2.act")
        test_data.act[[f]]$predict.g3.act <- predict(model.g3.act[[f]], newdata=test_data.act[[f]][,1:(2 + ncol(X) + Ncol_g234)])
        print("test_data.act[[f]]$predict.g3.act")
        test_data.act[[f]]$predict.g4.act <- predict(model.g4.act[[f]], newdata=test_data.act[[f]][,1:(2 + ncol(X) + Ncol_g234)])
        print("test_data.act[[f]]$predict.g4.act")
        test_data.act[[f]]$predict.m.act  <- predict(model.m.act[[f]],  newdata=test_data.act[[f]][,1:(2 + ncol(X))])
        print("test_data.act[[f]]$predict.m.act")
        mse_per_fold.g2.act[[f]] <- caret::postResample(pred = test_data.act[[f]]$predict.g2.act, obs = test_data.act[[f]]$Y)
        mse_per_fold.g3.act[[f]] <- caret::postResample(pred = test_data.act[[f]]$predict.g3.act, obs = test_data.act[[f]]$Y)
        mse_per_fold.g4.act[[f]] <- caret::postResample(pred = test_data.act[[f]]$predict.g4.act, obs = test_data.act[[f]]$Y)
        mse_per_fold.m.act[[f]]  <- caret::postResample(pred = test_data.act[[f]]$predict.m.act,  obs = test_data.act[[f]]$Y)
      } else if (length(unique(Y)) == 2 ){
        test_data.act[[f]]$predict.g2.act <-  as.numeric(predict(model.g2.act[[f]], newdata=test_data.act[[f]][,1:(2 + ncol(X) + Ncol_g234)], type="response"))
        print("test_data.act[[f]]$predict.g2.act")
        test_data.act[[f]]$predict.g3.act <-  as.numeric(predict(model.g3.act[[f]], newdata=test_data.act[[f]][,1:(2 + ncol(X) + Ncol_g234)], type="response"))
        print("test_data.act[[f]]$predict.g3.act")
        test_data.act[[f]]$predict.g4.act <-  as.numeric(predict(model.g4.act[[f]], newdata=test_data.act[[f]][,1:(2 + ncol(X) + Ncol_g234)], type="response"))
        print("test_data.act[[f]]$predict.g4.act")
        test_data.act[[f]]$predict.m.act  <-  as.numeric(predict(model.m.act[[f]],  newdata=test_data.act[[f]][,1:(2 + ncol(X))],             type="response"))
        print("test_data.act[[f]]$predict.m.act")
        test_data.act[[f]] <- test_data.act[[f]] %>% dplyr:: mutate(Y_predict.g2.act = factor(ifelse(predict.g2.act >0.5, 1, 0) ), 
                                                                    Y_predict.g3.act = factor(ifelse(predict.g3.act >0.5, 1, 0) ), 
                                                                    Y_predict.g4.act = factor(ifelse(predict.g4.act >0.5, 1, 0) ), 
                                                                    Y_predict.m.act  = factor(ifelse(predict.m.act >0.5, 1, 0) ), 
                                                                    Y=factor(Y))
        print("test_data.act[[f]]")
        ConfusionM_g2.act <- caret::confusionMatrix(data = test_data.act[[f]]$Y_predict.g2.act, reference = test_data.act[[f]]$Y, positive = "1")
        print(" ConfusionM_g2.act")
        ConfusionM_g3.act <- caret::confusionMatrix(data = test_data.act[[f]]$Y_predict.g3.act, reference = test_data.act[[f]]$Y, positive = "1")
        print(" ConfusionM_g3.act")
        ConfusionM_g4.act <- caret::confusionMatrix(data = test_data.act[[f]]$Y_predict.g4.act, reference = test_data.act[[f]]$Y, positive = "1")
        print(" ConfusionM_g4.act")
        ConfusionM_m.act  <- caret::confusionMatrix(data = test_data.act[[f]]$Y_predict.m.act,  reference = test_data.act[[f]]$Y, positive = "1")
        print(" ConfusionM_m.act")
        accuracy_per_fold.g2.act[[f]] <- as.numeric(ConfusionM_g2.act[[3]][1])
        print("accuracy_per_fold.g2.act[[f]]")
        accuracy_per_fold.g3.act[[f]] <- as.numeric(ConfusionM_g3.act[[3]][1])
        print("accuracy_per_fold.g3.act[[f]]")
        accuracy_per_fold.g4.act[[f]] <- as.numeric(ConfusionM_g4.act[[3]][1])
        print("accuracy_per_fold.g4.act[[f]]")
        accuracy_per_fold.m.act[[f]]  <- as.numeric(ConfusionM_m.act[[3]][1])
        print("accuracy_per_fold.m.act[[f]]")
      }
      #transformed outcome will always be continuous 
      test_data.tran[[f]]$predict.g2.tran <- predict(model.g2.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g234)])
      print("test_data.tran[[f]]$predict.g2.tran")
      test_data.tran[[f]]$predict.g3.tran <- predict(model.g3.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g234)])
      print("test_data.tran[[f]]$predict.g3.tran")
      test_data.tran[[f]]$predict.g4.tran <- predict(model.g4.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g234)])
      print("test_data.tran[[f]]$predict.g4.tran")
      test_data.tran[[f]]$predict.m.tran  <- predict(model.m.tran[[f]],  newdata=test_data.tran[[f]][,1:(2 + ncol(X))])
      print("test_data.tran[[f]]$predict.m.tran")
      mse_per_fold.g2.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g2.tran, obs = test_data.tran[[f]]$Y)
      print("mse_per_fold.g2.tran[[f]]")
      mse_per_fold.g3.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g3.tran, obs = test_data.tran[[f]]$Y)
      print("mse_per_fold.g3.tran[[f]]")
      mse_per_fold.g4.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g4.tran, obs = test_data.tran[[f]]$Y)
      print("mse_per_fold.g4.tran[[f]]")
      mse_per_fold.m.tran[[f]]  <- caret::postResample(pred = test_data.tran[[f]]$predict.m.tran,  obs = test_data.tran[[f]]$Y)
      print("mse_per_fold.m.tran[[f]]")
    } #f loop over
    ##############################  
    # f iteration ends
    ############################## 
    P_all_folds <- as.data.frame(do.call("rbind", HTE_P_cf.raw.L))
    print("P_all_folds")
    P_all_folds.mean <- mean(P_all_folds$V1)
    print("P_all_folds.mean")
    #transformed outcome will always be continuous 
    mse_all_folds.g2.tran <- CVBIAS_D_MAJORITY  (vote_D2_subgroup.L, mse_per_fold.g2.tran, "MSE")
    print(" mse_all_folds.g2.tran")
    mse_all_folds.g3.tran <- CVBIAS_D_MAJORITY  (vote_D3_subgroup.L, mse_per_fold.g3.tran, "MSE")
    print(" mse_all_folds.g3.tran")
    mse_all_folds.g4.tran <- CVBIAS_D_MAJORITY  (vote_D4_subgroup.L, mse_per_fold.g4.tran, "MSE")
    print(" mse_all_folds.g4.tran")
    mse_all_folds.m.tran  <- CVBIAS_MAIN  (mse_per_fold.m.tran, "MSE")
    print(" mse_all_folds.m.tran")
    grp_pick.tran=which.min(c(mse_all_folds.g2.tran, mse_all_folds.g3.tran, mse_all_folds.g4.tran))
    print("grp_pick.tran")
    #grp_CV_pick_D.tran  <- PICK_m234(P_all_folds.mean, mse_all_folds.m.tran, mse_all_folds.g2.tran, mse_all_folds.g3.tran, mse_all_folds.g4.tran)
    
    # grp_pick.tran is already based on majority voted value (if all are unique, then pick the subgroup decision with smallest group number)
    #  CV_SG_MAJORITY will choose subgroup decision again following 3 scenarios in the CVBIAS_D_MAJORITY function
    seletedSG_CV.tran.bias <- CV_SG_MAJORITY(  eval(parse(text =  paste0("vote_D",  (grp_pick.tran + 1), "_subgroup.L" ))) ) 
    print("seletedSG_CV.tran.bias")
    if (length(unique(Y)) >=8){
      # Bind together, add MSE
      mse_all_folds.m.act  <- CVBIAS_MAIN  (mse_per_fold.m.act, "MSE")
      print("mse_all_folds.m.act")
      mse_all_folds.g2.act  <- CVBIAS_D_MAJORITY  (vote_D2_subgroup.L, mse_per_fold.g2.act, "MSE")
      print("mse_all_folds.g2.act")
      mse_all_folds.g3.act  <- CVBIAS_D_MAJORITY  (vote_D3_subgroup.L, mse_per_fold.g3.act, "MSE")
      print("mse_all_folds.g3.act")
      mse_all_folds.g4.act  <- CVBIAS_D_MAJORITY  (vote_D4_subgroup.L, mse_per_fold.g4.act, "MSE")
      print("mse_all_folds.g4.act")
      grp_pick.act =which.min(c( mse_all_folds.g2.act,  mse_all_folds.g3.act,  mse_all_folds.g4.act))
      print("grp_pick.act") #grp_CV_pick_D.act   <- PICK_m234(P_all_folds.mean, mse_all_folds.m.act,  mse_all_folds.g2.act,  mse_all_folds.g3.act,  mse_all_folds.g4.act)
    } else if (length(unique(Y)) ==2){
      accuracy_all_folds.m.act   <- CVBIAS_MAIN  (accuracy_per_fold.m.act, "accuracy")
      print("accuracy_all_folds.m.act") 
      accuracy_all_folds.g2.act  <- CVBIAS_D_MAJORITY  (vote_D2_subgroup.L, accuracy_per_fold.g2.act, "accuracy")
      print("accuracy_all_folds.g2.act") 
      accuracy_all_folds.g3.act  <- CVBIAS_D_MAJORITY  (vote_D3_subgroup.L, accuracy_per_fold.g3.act, "accuracy")
      print("accuracy_all_folds.g3.act") 
      accuracy_all_folds.g4.act  <- CVBIAS_D_MAJORITY  (vote_D4_subgroup.L, accuracy_per_fold.g4.act, "accuracy")
      print("accuracy_all_folds.g4.act") 
      grp_pick.act =which.min(c(accuracy_all_folds.g2.act,  accuracy_all_folds.g3.act,  accuracy_all_folds.g4.act))
      print("grp_pick.act") 
      #grp_CV_pick_D.act   <- PICK_m234(P_all_folds.mean, accuracy_all_folds.m.act,  accuracy_all_folds.g2.act,  accuracy_all_folds.g3.act, accuracy_all_folds.g4.act)
    }
    seletedSG_CV.act.bias  <- CV_SG_MAJORITY(  eval(parse(text =  paste0("vote_D",  (grp_pick.act  + 1),  "_subgroup.L" ))) )
    print("seletedSG_CV.act.bias") 
    #if majority (IFM)
    IFM_Deci_AIC_test.act  <-lapply( Deci_AIC_test.act,  
                                     function(df) 
                                    identical(df, unique(Deci_AIC_test.act)[[which.max(MAJORITY_COUNT(unique(Deci_AIC_test.act),Deci_AIC_test.act))]]  )) 
    print("IFM_Deci_AIC_test.act") 
    IFM_Deci_AIC_test.tran <-lapply( Deci_AIC_test.tran, 
                                     function(df) 
                                    identical(df, unique(Deci_AIC_test.tran)[[which.max(MAJORITY_COUNT(unique(Deci_AIC_test.tran),Deci_AIC_test.tran))]] ))  
    print("IFM_Deci_AIC_test.tran") 
    seletedSG_CV.act.AIC  <- STANDARD_CHECK(IFM_Deci_AIC_test.act,  Deci_AIC_test.act)[[1]]
    print("seletedSG_CV.act.AIC") 
    seletedSG_CV.tran.AIC <- STANDARD_CHECK(IFM_Deci_AIC_test.tran, Deci_AIC_test.tran)[[1]]
    print("seletedSG_CV.tran.AIC") 
    seletedSG_NCV.act.AIC <- Deci_AIC_test.act[[1]]
    print("seletedSG_NCV.act.AIC") 
    seletedSG_NCV.tran.AIC <- Deci_AIC_test.tran[[1]]
    print("seletedSG_NCV.AIC.AIC") 
  } else if (round(HTE_P_cf.raw,1) > SigLevel){
    CVL_Deci_AIC_test.act="NA"
    CVL_Deci_AIC_test.tran="NA"
    seletedSG_CV.act.bias  = "NA"
    seletedSG_CV.tran.bias = "NA"
    seletedSG_CV.act.AIC  = "NA"
    seletedSG_CV.tran.AIC = "NA"
    seletedSG_NCV.act.AIC  = "NA"
    seletedSG_NCV.tran.AIC = "NA"
    vote_D4_subgroup.L="NA"
    vote_D3_subgroup.L="NA"
    vote_D2_subgroup.L="NA"
  }
  e_iCF_CV=Sys.time()
  time_iCF_CV = difftime(e_iCF_CV, s_iCF_CV, units="secs")
  return(list( CVL_Deci_AIC_test.act =  Deci_AIC_test.act,
               CVL_Deci_AIC_test.tran = Deci_AIC_test.tran,
              seletedSG_CV.act.bias  = seletedSG_CV.act.bias ,
              seletedSG_CV.tran.bias = seletedSG_CV.tran.bias,
              seletedSG_CV.act.AIC  = seletedSG_CV.act.AIC ,
              seletedSG_CV.tran.AIC = seletedSG_CV.tran.AIC,
              seletedSG_NCV.act.AIC  = seletedSG_NCV.act.AIC ,
              seletedSG_NCV.tran.AIC = seletedSG_NCV.tran.AIC,
              vote_D4_subgroup.L= vote_D4_subgroup.L,
              vote_D3_subgroup.L= vote_D3_subgroup.L,
              vote_D2_subgroup.L= vote_D2_subgroup.L,
              leafsize =  leafsize,
              time_iCFCV = time_iCF_CV))
  
}





