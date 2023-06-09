####################################################################
# SUBGROUP DECISION 
# Author: Tiansheng Wang
# Last update date: 12/8/2020
# Version: 0.1      
####################################################################

#' Wrapper function that convert tree structure to subgroup 
#' @param tree tree structure df
#'  
#' @return subgroup in list format
#' 
#' @export

TREE2SUBGROUP <- function(tree){
  subgroup.ori    <- PRE_MAJORITY_SUBGROUP ( list(tree ) )
  
  subgroup_key.ori   <-lapply(subgroup.ori,   function(df) subset(df, select=c("subgroupID", "subgroup")))
  
  #did not vote for subgroup, use the synthetic voted tree structure to get subgroup decision
  vote_subgroup.ori <- list(majority=subgroup_key.ori[[1]])
  
  return(vote_subgroup.ori )
}


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


#' Function that does prepare data that
#' 1) ready to apply the formulas including interacting terms between treatment and factor subgroup deicision variable; 
#' 2)transformed outcome
#' @param dat P-value for raw full CF
#' @param outcome_type "actual" or "transform" 
#'  
#' @return the prepared data and "contr"
#' 
#' @export
#' https://gsbdbi.github.io/ml_tutorial/hte_tutorial/hte_tutorial.html
#' Estimation of Heterogeneous Treatment Effects - prepared for “Machine Learning and Causal Inference” class

SGMODEL_DATA<-function(dat, outcome_type_at){
  if (intTRUE != "Unknown" & intTRUE != "Medicare"){
    #remove the columns for subgroup definitions (characters)
    dat_ID_SG_pre <- dat %>%  dplyr::select(-contains( c("G5_define", "G4_define", "G3_define", "G2_define" ,"ID" ) ) ) 
  } else {
    #for Medicare, don't use contains Fx as it delete more columns than expected!!!
    dat_ID_SG_pre <- dat %>%  # dplyr::select(-contains( c("G4_define", "G3_define", "G2_define" ,"ID" ) ) ) #%>% 
      dplyr::select(-G5_define, -G4_define, -G3_define, -G2_define ,-ID  ) 
  } 
  
  dat_ID_SG_pre$ps <- prop.func(as.matrix( dat_ID_SG_pre %>% dplyr::select (-contains( c("W", "Y", "G5",  "G4", "G3", "G2")))), #matrix of covariates
                                dat_ID_SG_pre$W #outcome for glmnet
                                )
  #create transformed outcome, https://johaupt.github.io/blog/Uplift_ITE_summary.html
  dat_ID_SG_pre$Y_star <- ifelse(dat_ID_SG_pre$W==1, dat_ID_SG_pre$Y/dat_ID_SG_pre$ps, -dat_ID_SG_pre$Y/(1-dat_ID_SG_pre$ps))
  
  # Compute Y-star, https://gsbdbi.github.io/ml_tutorial/hte_tutorial/hte_tutorial.html 
  #dat_ID_SG_pre$p <- mean(dat_ID_SG_pre$W)
  #dat_ID_SG_pre$Y_star <- ((dat_ID_SG_pre$W - dat_ID_SG_pre$p)/(dat_ID_SG_pre$p*(1-dat_ID_SG_pre$p)))*dat_ID_SG_pre$Y
  

  #create dataset depends on mehtod (actual outcome vs transormed outcome)
  if (outcome_type_at=="actual"){
    dat_ID_SG_pre <- dat_ID_SG_pre %>% dplyr::select(-c(ps, Y_star))
  } else if (outcome_type_at=="transform"){
    
    dat_ID_SG_pre <- dat_ID_SG_pre %>% dplyr::select(-c(ps, Y)) %>%
      dplyr::rename(Y=Y_star) %>% #rename Y_start as Y
      dplyr::select(Y, everything())
  }
  
  
  if (identical(vars_catover2,NA) ){
    dat_ID_SG_df <-  dat_ID_SG_pre %>%  mutate_at(vars( c ("G5", "G4", "G3", "G2")    ), as.factor) 
    contr <- rep(list("contr.sum"), ncol( dat_ID_SG_df %>% dplyr::select( c("G5", "G4","G3", "G2") )  
                                         )
                                         )
    names(contr) <-  names( dat_ID_SG_df %>% dplyr::select(c("G5", "G4", "G3", "G2"  )) 
                            )   # grep("_", names(dat_ID_SG %>% select( -contains(c("_define")))), value=TRUE)
  } else {
    dat_ID_SG_df <-  dat_ID_SG_pre %>%  mutate_at(vars(  c ("G5", "G4", "G3", "G2", vars_catover2) ) , as.factor) 
    contr <- rep(list("contr.sum"), ncol( dat_ID_SG_df %>% dplyr::select(   c("G5", "G4","G3", "G2",vars_catover2 ) ) ) )
    names(contr) <-  names( dat_ID_SG_df %>% dplyr::select(c("G5", "G4", "G3", "G2", vars_catover2  ) ) )   # grep("_", names(dat_ID_SG %>% select( -contains(c("_define")))), value=TRUE)
  } 
  
  return(list(dat_ID_SG_df=dat_ID_SG_df,
              contr=contr))
  
}

#' Function that does shows the depth of subgroup selected 
#' @param HTE_P_value P-value for raw full CF
#' @param m   metric for main effects model 
#' @param g2  metric for g2 model   
#' @param g3  metric for g3 model 
#' @param g4  metric for g4 model 
#' @param g5  metric for g5 model 
#'  
#' @return the depth of subgroup selected, e.g. "W:G2x"
#' 
#' @export

PICK_m2345 <- function(HTE_P_value, m, g2, g3, g4, g5, P_threshold){
  if (round(HTE_P_value,1) >  P_threshold){#may be false negative, i.e. HTE may exist, 
    grp_pick =  paste0("W:G",  which.min(c(m, g2, g3, g4, g5)), "x" )
    grp_pick_D = ifelse(grp_pick=="W:G1x", "NA",  grp_pick)
  } else {#if HTE P-value <0.1, then we know for sure HTE exist, 
    #i.e. to reduce subgroup number and number of covariates define a subgroup
    grp_pick_D= paste0("W:G",  which.min(c(g2, g3, g4, g5)) + 1 , "x" )
  }
  return(grp_pick_D)
}


#' Function that does shows the depth of subgroup selected 
#' @param HTE_P_value P-value for raw full CF
#' @param MAF_m   the mean MSE value across CV for main effect model 
#' @param MAF_g2  the mean MSE value across CV for D2 SG model   
#' @param MAF_g3  the mean MSE value across CV for D3 SG model 
#' @param MAF_g4  the mean MSE value across CV for D4 SG model 
#' @param MAF_g5  the mean MSE value across CV for D5 SG model
#' @param V_D2_SG.L  the mean MSE value across CV for D2 SG model   
#' @param V_D3_SG.L  the mean MSE value across CV for D3 SG model 
#' @param V_D4_SG.L  the mean MSE value across CV for D4 SG model 
#' @param V_D5_SG.L  the mean MSE value across CV for D5 SG model
#' @param S_D2_T_r  the mean MSE value across CV for D2 SG model   
#' @param S_D3_T_r  the mean MSE value across CV for D3 SG model 
#' @param S_D4_T_r  the mean MSE value across CV for D4 SG model 
#' @param S_D5_T_r  the mean MSE value across CV for D5 SG model
#' @param S_D2_CV  the mean MSE value across CV for D2 SG model   
#' @param S_D3_CV  the mean MSE value across CV for D3 SG model 
#' @param S_D4_CV  the mean MSE value across CV for D4 SG model 
#' @param S_D5_CV  the mean MSE value across CV for D5 SG model
#' @param S_D2_vote  the mean MSE value across CV for D2 SG model   
#' @param S_D3_vote  the mean MSE value across CV for D3 SG model 
#' @param S_D4_vote  the mean MSE value across CV for D4 SG model 
#' @param S_D5_vote  the mean MSE value across CV for D5 SG model
#' @param P_threshold  P-value threshold 
#' @param K  cross-validation fold 
#'  
#' @return selected SG and it's stability of CV, and stability of Vote
#' 
#' @export

PICK_CV <- function(HTE_P_value, MAF_m, 
                    MAF_g2, MAF_g3, MAF_g4, MAF_g5, 
                    V_D2_SG.L, V_D3_SG.L, V_D4_SG.L, V_D5_SG.L,
                    S_D2_T_r, S_D3_T_r, S_D4_T_r, S_D5_T_r,
                    #S_D2_CV, S_D3_CV, S_D4_CV, S_D5, CV,
                    #S_D2_vote, S_D3_vote, S_D4_vote, S_D5, vote,
                    P_threshold, K){

  grp_pick        = which.min(c(MAF_g2,       MAF_g3,        MAF_g4,        MAF_g5))

  
  #Function that does shows the depth of subgroup selected 
  #grp_CV_pick_D.act   <- PICK_m2345(HTE_P_cf.raw, #P_all_folds.mean, 
  #                                  mse_all_folds.m.act,  mse_all_folds.g2.act,  mse_all_folds.g3.act,  mse_all_folds.g4.act, mse_all_folds.g5.act, P_threshold)
  grp_pick_D <- PICK_m2345(HTE_P_cf.raw, MAF_m, MAF_g2, MAF_g3, MAF_g4, MAF_g5, P_threshold)
  
  
  #finds the majority voted subgroup decision across CV
  #seletedSG_CV.act.bias  <- CV_SG_MAJORITY(  eval(parse(text =  paste0("vote_D",  grp_pick.act  + 1, "_subgroup.L" ))) )
  selectedSG <- CV_SG_MAJORITY( eval(parse(text =  paste0("vote_D",      grp_pick + 1, "_subgroup.L" ))),
                                eval(parse(text =  paste0("stability_D", grp_pick + 1, "_T_r" ))),
                                K)$majority_CV 
  
  selectedSG_stability_CV <- #CV_SG_MAJORITY(eval(parse(text =  paste0("vote_D",      grp_pick + 1, "_subgroup.L" ))),
                                    #              eval(parse(text =  paste0("stability_D", grp_pick + 1, "_T_r" ))),
                                    #              K)$stability_CV 
                              eval(parse(text =  paste0("stability_D", grp_pick + 1, "_CV" )))
  
  selectedSG_stability_vote <-  eval(parse(text =  paste0("stability_D", grp_pick + 1, "_vote" )))
  
  
  
  selectedSG_wt <- as.data.frame(cbind(grp_pick_D, selectedSG_stability_vote, selectedSG_stability_CV, 
                                             round(MAF_g2,5), 
                                             round(MAF_g3,5), 
                                             round(MAF_g4,5), 
                                             round(MAF_g5,5) )) %>%
                    dplyr::rename(grp_pick =grp_pick_D,
                                  stability_vote_mean4CV=selectedSG_stability_vote, 
                                  stability_CV= selectedSG_stability_CV,
                                  mse_g2_Model= V4,
                                  mse_g3_Model= V5,
                                  mse_g4_Model= V6,
                                  mse_g5_Model= V7) %>%
                    
                    knitr::kable()
  
            
  return(list( selectedSG = selectedSG ,
               selectedSG_wt   =   selectedSG_wt))
}


CATE_SG <- function(SGdecision){
  Delta_Y <- DltY_DATA(SGdecision, Train_ID, "non-test")
  
  CATE_table <- Delta_Y %>% dplyr::select(SubgroupID, Definition, W, Y,
                                          dltY_crude, dltY_crude_low, dltY_crude_up, dltY_ate,dltY_ate_low,dltY_ate_up, dltY_att, dltY_att_low, dltY_att_up) %>%
    dplyr::group_by(SubgroupID, Definition,W,Y) %>%
    dplyr::summarise(n= n(),
                     dltY_crude = mean(dltY_crude),
                     dltY_crude_low = mean(dltY_crude_low),
                     dltY_crude_up = mean(dltY_crude_up),
                     dltY_ate = mean(dltY_ate),
                     dltY_ate_low = mean(dltY_ate_low),
                     dltY_ate_up = mean(dltY_ate_up),
                     dltY_att = mean(dltY_att),
                     dltY_att_low= mean(dltY_att_low),
                     dltY_att_up = mean(dltY_att_up))%>%
    dplyr::arrange(dltY_ate)%>%
    dplyr::mutate(#Treatment = sum(W)
      CATE_crude=paste0(dltY_crude,"(",dltY_crude_low,",",dltY_crude_up,")"),
      CATE_iptw  =paste0(dltY_ate,"(",dltY_ate_low,",",dltY_ate_up,")"),
      CATE_smr  =paste0(dltY_att,"(",dltY_att_low,",",dltY_att_up,")")) %>%
    dplyr::select(-c(dltY_crude, dltY_crude_low, dltY_crude_up, dltY_ate, dltY_ate_low, dltY_ate_up, dltY_att, dltY_att_low, dltY_att_up)) 
  
  return(CATE_table)
  
}





#' Function that does shows the final subgroup decisions 
#' @param grp_pick P-value for raw full CF
#' @param iterationNo   
#'  
#' @return the final subgroup decision,
#' 
#' @export
FINALDECI<- function(grp_pick, V_D4_subgroup,  V_D3_subgroup, V_D2_subgroup){
  #if (iterationNo > 1){
  CFmethod_prefix  = "vote_D"
  CFmethod_postfix = "_subgroup"
  #} else if (iterationNo ==1 ){
  #  CFmethod_prefix  =  "vote_D" #"vote_1D"
  #  CFmethod_postfix = "_subgroup" #"b_subgroup"
  #
  V_list=list(D2= V_D2_subgroup,
              D3= V_D3_subgroup,
              D4= V_D4_subgroup)
  if (grp_pick != "NA"){
    #Deci_Final        <- eval(parse(text =  paste0(CFmethod_prefix,  stringr::str_sub(grp_pick,     4,4), CFmethod_postfix )))$majority
    Deci_Final        <-  V_list[[as.numeric( stringr::str_sub(grp_pick,     4,4)) - 1]]$majority
  } else {
    Deci_Final = "NA"
  }
  return(Deci_Final)
}



#' #Function that obtain MODELS from subgroup decisions at different depths from CF    
#' @param dat the dataset for group lasso 
#' @param HE_P_value raw heterogeneity test P-value 
#' @param formula_gl formula for group lasso 
#' @param method_CF method for CF: "iCF" or "oneCF"
#' @param method_CF actall outcome or transformed outcome, iCFCV funciton will proivde results for both
#' @param method_deci method for selecting "cv" or "shrink"  
#'  
#' @return the data for group lasso
#' 
#' @export


CF_GROUP_DECISION <- function(HTE_P_cf.raw, dat, method_CF, outcome_type_at ,  P_threshold ){
  if (  round(HTE_P_cf.raw,1) > P_threshold ){
    model.basic="NA"
    model.g2="NA"
    model.g3="NA"
    model.g4="NA"
    model.g5="NA"
    return(list(     
                  model.basic=model.basic,
                  model.g2=model.g2,
                  model.g3=model.g3,
                  model.g4=model.g4,
                  model.g5=model.g5))
  } else {
    #################################f
    ##      PREPARATION            ##
    #################################
    #Function that does prepare data that 
    #1) ready to apply the formulas including interacting terms between treatment and factor subgroup deicision variable; 
    #2)transformed outcome
    SGmodel_dat <-  SGMODEL_DATA (dat,outcome_type_at)#works for HD covariates!
    
    dat_ID_SG_df = SGmodel_dat$dat_ID_SG_df
    contr = SGmodel_dat$contr
    
    
    if (length(unique(dat_ID_SG_df$Y)) >=8 #transformed outcome
        ){
      model.basic = lm(formula_basic, data = dat_ID_SG_df) 
      if ("G2" %in% names(contr)) {model.g2    = lm(formula_g2,    data = dat_ID_SG_df)}
      if ("G3" %in% names(contr)) {model.g3    = lm(formula_g3,    data = dat_ID_SG_df)}
      if ("G4" %in% names(contr)) {model.g4    = lm(formula_g4,    data = dat_ID_SG_df)}
      if ("G5" %in% names(contr)) {model.g5    = lm(formula_g5,    data = dat_ID_SG_df)}
    } else if (length(unique(dat_ID_SG_df$Y)) ==2) { #actual binary outcome,
      #binary outcome
      model.basic = glm(formula_basic, data = dat_ID_SG_df, family=binomial)
      if ("G2" %in% names(contr)) {model.g2    = glm(formula_g2,    data = dat_ID_SG_df, family=binomial)}
      if ("G3" %in% names(contr)) {model.g3    = glm(formula_g3,    data = dat_ID_SG_df, family=binomial)}
      if ("G4" %in% names(contr)) {model.g4    = glm(formula_g4,    data = dat_ID_SG_df, family=binomial)}
      if ("G5" %in% names(contr)) {model.g5    = glm(formula_g5,    data = dat_ID_SG_df, family=binomial)}
    } 
    
    return(list(
      model.basic=model.basic,
      model.g2=model.g2,
      model.g3=model.g3,
      model.g4=model.g4,
      model.g5=model.g5))
  } #use HTE p-value at the begnining
}






#------------------------------------------------------------------------------------------------------------
# IV. ACCURACY (=TRUE OR NOT) OF EACH TYPE OF TREE AND FINAL SUBGROUP DECISION BY DIFFERENT THRESHOLD VALUE 
#------------------------------------------------------------------------------------------------------------
TEST_DECISION <- function(test_df, true_df){
  if (identical(test_df, true_df)){ 
    test.score <- 1
  } else {
    test.score <- 0
  }
  return(test.score)
}

#------------------------------------------------------------------------------------------------------------
# IV. Which Depth iCF decision comes from 
#------------------------------------------------------------------------------------------------------------
WHICH_DECISION <- function(depth_df, final_df){
  if (identical(depth_df, final_df)){ 
    test.score <- 1
  } else {
    test.score <- 0
  }
  return(test.score)
}

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
TorF <- function(v_D4_tree, v_D4_tree_R, Desc_v_D4, v_D3_tree, v_D3_tree_R, Desc_v_D3, v_D2_tree, v_D2_tree_R, Desc_v_D2, truth.list){
  TF_D4_T     <- ifelse (identical(v_D4_tree$majority,      truth.list$tree_true), 1, 0)
  TF_D4_T_N1  <- ifelse (identical(v_D4_tree$majority[1,],  truth.list$tree_true_N1), 1, 0)   #get from tree list, cannot check from subgroup or tree_r list
  TF_D4_T_N123<- ifelse (identical(v_D4_tree$majority[1:3,],truth.list$tree_true_N123), 1, 0) #get from tree list, cannot check from subgroup or tree_r list
  TF_D4_T_r   <- ifelse (identical(v_D4_tree_R$majority,    truth.list$tree_true_r), 1, 0)
  TF_D4_SG    <- ifelse (identical(Desc_v_D4,               truth.list$tree_true_subgroup), 1, 0)
  
  TF_D3_T     <- ifelse (identical(v_D3_tree$majority,      truth.list$tree_true), 1, 0)
  TF_D3_T_N1  <- ifelse (identical(v_D3_tree$majority[1,],  truth.list$tree_true_N1), 1, 0)   #get from tree list, cannot check from subgroup or tree_r list
  TF_D3_T_N123<- ifelse (identical(v_D3_tree$majority[1:3,],truth.list$tree_true_N123), 1, 0) #get from tree list, cannot check from subgroup or tree_r list
  TF_D3_T_r   <- ifelse (identical(v_D3_tree_R$majority,    truth.list$tree_true_r), 1, 0)
  TF_D3_SG    <- ifelse (identical(Desc_v_D3,               truth.list$tree_true_subgroup), 1, 0)
  
  TF_D2_T     <- ifelse (identical(v_D2_tree$majority,      truth.list$tree_true), 1, 0)
  TF_D2_T_N1  <- ifelse (identical(v_D2_tree$majority[1,],  truth.list$tree_true_N1), 1, 0)   #get from tree list, cannot check from subgroup or tree_r list
  TF_D2_T_N123<- ifelse (identical(v_D2_tree$majority[1:3,],truth.list$tree_true_N123), 1, 0) #always wrong (i.e 0) when truth is "deeper tree"
  TF_D2_T_r   <- ifelse (identical(v_D2_tree_R$majority,    truth.list$tree_true_r), 1, 0)
  TF_D2_SG    <- ifelse (identical(Desc_v_D2,               truth.list$tree_true_subgroup), 1, 0)
  
  return(
    list(  TF_D4_T     = TF_D4_T, 
           TF_D4_T_N1  = TF_D4_T_N1,
           TF_D4_T_N123= TF_D4_T_N123,
           TF_D4_T_r   = TF_D4_T_r ,
           TF_D4_SG    = TF_D4_SG ,
           
           TF_D3_T     = TF_D3_T,
           TF_D3_T_N1  = TF_D3_T_N1,
           TF_D3_T_N123= TF_D3_T_N123, 
           TF_D3_T_r   = TF_D3_T_r ,
           TF_D3_SG    = TF_D3_SG,
           
           TF_D2_T     = TF_D2_T, 
           TF_D2_T_N1  = TF_D2_T_N1 ,
           TF_D2_T_N123= TF_D2_T_N123 , 
           TF_D2_T_r   = TF_D2_T_r ,
           TF_D2_SG    = TF_D2_SG  
    )
  )
}


N_SUBGROUP <- function(deci, method_L) {
  
  if (is.null (nrow(deci))==TRUE  & method_L == "tree_sg" ) {
    n_subgroup = 0
  }  else if (is.null (nrow(deci))==FALSE & method_L == "tree_sg") {
    n_subgroup = nrow(deci)
  }  else if (deci[1] =="NA"              & method_L == "non-tree_int" ){
    n_subgroup = 0
  }  else if (deci[1] !="NA"              & method_L == "non-tree_int" ){
    n_subgroup = length(deci)
  }
  return(n_subgroup)
}


DISCOVERY <- function(N_method, N_truth, TF){
  discovery="NA"
  if (N_method < N_truth){
    discovery="typeII"
  } else if (N_method > N_truth ) {
    discovery="typeI"
  } else if (N_method == N_truth & TF==1) {
    discovery="accurate"
  }
  return(discovery)
}


#https://community.rstudio.com/t/detect-presence-absence-of-elements-in-nested-lists/31515
#' Function check if vector string in 2nd position nested in the list in 1st position
#' @param your_list list
#' @param vector_strings vector of string
#' @return T of F
#' 
#' @export
#' 
DETECT_STRING <- function(your_list, vector_strings){
  lapply(your_list, function(x) {
    if(TRUE %in% str_detect(x, paste(vector_strings, collapse = "|"))){
      TRUE
    } else {FALSE}
  })
}

#' Function that assess accuracy, false positive, and false negative.
#' Accuracy definition: all subgup decisions (for tree-based methods that explicitly give subgroup) or all interactions (for non-tree based methds) are accurately identified  
#' @param decision the final decision could be subgroup dataframe or interction vector
#' @param truth_description the description of truth
#' @param tree_true_subgroup the true subgroup
#' @param truth_INT the true interaction, when the metric==CATEmax, use true CATEmax subgroup instead of true interaction!!!
#' @param method different methods (tree-based or non-tree-based)
#' @param metric assesing perforamnce of methods by subgroup or interaction identification
#' 
#' @return value=score of accuracy assessment; sg_Nacc_Nall=% of true subgroups among identified subgroups; sg_Nacc_NallT=% of true subgroups among total true subgroups 
#' 
#' @export
#' 

#disco_SG_AIC_iCF.ori      <- DISCOVERYRATE(Deci_AIC_iCF.ori,      truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")
#decision=Deci_AIC_iCF.ori
#method="iCF"
#metric="subgroup"

DISCOVERYRATE <- function (decision, truth_description, tree_true_subgroup, truth_INT, method, metric){
  
  #if assessing performance by subgroup identification for tree-based method:
  if ( method %in% c("oneCFv", "oneCFb", "iCF", "iCFv", "IT", "VT")==T & metric=="subgroup"){
    if (is.null( dim(decision))== F ) {
      Deci_vec    <- dplyr::pull(decision,           subgroup)#subgroup here is a column's name
    } else {
      Deci_vec = "NA"
    }
    
    if (is.null( dim(tree_true_subgroup)) == F) {
      Truth_vec  <- dplyr::pull(tree_true_subgroup, subgroup)
    } else {
      Truth_vec = "NA"
    }
    #get N to judge false+ or false-
    N_true   <- N_SUBGROUP (tree_true_subgroup, "tree_sg")
    N_deci   <- N_SUBGROUP (decision,           "tree_sg")
  } else if ( metric=="interaction" | metric=="CATEmax") {
    #if assessing performance by interaction identification for ALL method:
    
    
    Deci_vec       <- decision
    Truth_vec      <- truth_INT
    N_true  <- N_SUBGROUP (truth_INT,      "non-tree_int")
    N_deci  <- N_SUBGROUP (decision,       "non-tree_int")
  }
  
  Deci_vec_L  <- as.list(Deci_vec)  #the same for both TREE and non-TREE method 
  Truth_vec_L <- as.list(Truth_vec) #the same for both TREE and non-TREE method
  
  value = ifelse (
    ###########################--------------------------       
    # I.1 ACCURATE
    ###########################--------------------------
    (
      #----------------------------------------------------                                       
      #I.1.1. identical for tree-based method
      #----------------------------------------------------  
      ( metric=="subgroup"                           & identical(decision, tree_true_subgroup)  ) |
        ((metric=="CATEmax" |  metric=="interaction" ) & identical(decision, truth_INT #when metric=="CATEmax", truth_INT parameter is true subgroup with max CATE rather than overall true interaction
        )  )  |
        #----------------------------------------------------                                       
      #I.1.2. 4-way INT
      #----------------------------------------------------
      (metric=="interaction" & truth_description == "4-way INT of W, X3, X1, X8" & 
         #I.1.2.1. the same ignoring order
         ( setequal(Deci_vec, truth_INT) == T  | 
             #I.1.2.2. as long as W:X1:X3:X8 is identified:
             ( "W:X1:X3:X8" %in% decision  & 
                 #all inconsistent interactions are one of the following:
                 all( decision[is.na(match(decision, truth_INT))] %in% c("W:X1:X3", "W:X1:X8", "W:X3:X8", "W:X1", "W:X3", "W:X8") )  
             )
         )
      )
      |
        (metric=="interaction" & (truth_description == "4-way INT of W, X3, X1, X2(4l)" |  #ordinal splitter 4 levels
                                    truth_description == "4-way INT of W, X3, X1, X2(3l)" |  #ordinal splitter 3 levels
                                    truth_description == "4-way INT of W, X3, X1, X2(l)") &  #cotinuous splitter
           #I.1.3.1. the same ignoring order
           ( setequal(Deci_vec, truth_INT) == T | 
               #I.1.3.2. as long as W:X1:X3 is identified:
               ("W:X1:X2:X3" %in% decision  & 
                  #all inconsistent interactions are one of the following:
                  all( decision[is.na(match(decision, truth_INT))] %in% c("W:X1:X3", "W:X1:X2", "W:X2:X3", "W:X1", "W:X3", "W:X2") )  
               ) 
           )
        )
      |
        #----------------------------------------------------                                       
      #I.1.3. 3-way INT
      #----------------------------------------------------
      (metric=="interaction" & truth_description == "3-way INT of W, X3, X1"     & #binary splitter
         #I.1.3.1. the same ignoring order
         ( setequal(Deci_vec, truth_INT) == T | 
             #I.1.3.2. as long as W:X1:X3 is identified:
             ("W:X1:X3" %in% decision  & 
                #all inconsistent interactions are one of the following:
                all( decision[is.na(match(decision, truth_INT))] %in% c("W:X1", "W:X3") )
             ) 
         )
      )
      |
        (metric=="interaction" & (truth_description == "3-way INT of W, X3, X2(4l)" |  #ordinal splitter 4 levels
                                    truth_description == "3-way INT of W, X3, X2(3l)" |  #ordinal splitter 3 levels
                                    truth_description == "3-way INT of W, X3, X2(l)") &  #cotinuous splitter
           #I.1.3.1. the same ignoring order
           ( setequal(Deci_vec, truth_INT) == T | 
               #I.1.3.2. as long as W:X1:X3 is identified:
               ("W:X2:X3" %in% decision  & 
                  #all inconsistent interactions are one of the following:
                  all( decision[is.na(match(decision, truth_INT))] %in% c("W:X2", "W:X3") )
               ) 
           )
        )
      |
        #----------------------------------------------------                                       
      #I.1.4. 2-way INT, two or three 2 way INTs
      #----------------------------------------------------
      (metric=="interaction" & (substr(truth_description,1,5) == "2-way" | #2-way INT
                                  truth_description == "Three 2-way INTs of W, X3, X1, X8" | #three 2-way INT
                                  truth_description == "Two 2-way INTs of W, X3, X1") #two 2-way INT
       & setequal(Deci_vec, truth_INT) ) 
    ), 
    1 
    , 
    ############################--------------------------       
    # I.2 "FALSE NEGATIVE" abuse term
    ############################--------------------------
    #SG/INT No. < truth     & at 1st glance, none of defined subgroups are the same as the truth;   & ALL SG/INT in the decision vector are nested in truth list, e.g. X3:X1 in X3:X1:X8   
    ifelse( N_deci < N_true  & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                        & all( unlist(DETECT_STRING(Truth_vec_L, Deci_vec)) ) == T, 2,  
            #SG/INT No. < truth     & at 1st glance, none of defined subgroups are the same as the truth;   & some or all of SG/INT in the decision vector are nested in truth list, e.g. X3:X1 in X3:X1:X8   
            ifelse( N_deci < N_true  & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                        & all( unlist(DETECT_STRING(Truth_vec_L, Deci_vec)) ) == F, 3,  
                    # subgroup No. < truth  & some of defined subgroups are the same as the true subgroups;         & some SG/INT in the truth list not included in the decision vector    
                    ifelse( N_deci < N_true  & length( Deci_vec[Deci_vec %in% Truth_vec]) > 0                        & (length( Truth_vec[Truth_vec  %in% Deci_vec]) < length(Truth_vec) ) , 4, 
                            ############################--------------------------       
                            # I.3 "FALSE POSITIVE" abuse term
                            ############################--------------------------       
                            # SG/INT No. > truth  & at 1st glance, none of SG/INT are the same as the truths;               & NONE of SG/INT in the truth vector are nested in decision list  (totally wrong)
                            ifelse( N_deci > N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == F, 5,  
                                    #SG/INT No. >= truth  & at 1st glance, none of SG/INT are the same as the truths;                & some or all of SG/INT in the truth vector (eg,"W:X1" "W:X3") are nested in decision list (eg, "W:X1:X10"    "W:X1:X10:X4" "W:X1:X10:X7"), e.g.
                                    ifelse( N_deci > N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == T, 6,  
                                            # SG/INT No. >= truth  & some SG/INT are the same as in the true subgroups;                      & some SG/INT in the decision list not included in the truth vector                         
                                            ifelse( N_deci > N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) > 0                        & ( length( Deci_vec[Deci_vec %in% Truth_vec]) < length(Deci_vec)  ) , 7, 
                                                    ############################--------------------------       
                                                    # I.4 
                                                    ############################-------------------------- 
                                                    # SG/INT No. = truth  & at 1st glance, none of SG/INT are the same as the truths;               & NONE of SG/INT in the truth vector are nested in decision list (totally wrong)
                                                    ifelse( N_deci == N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == F, 8,  
                                                            #SG/INT No. >= truth  & at 1st glance, none of SG/INT are the same as the truths;                & some or all of SG/INT in the truth vector (eg,"W:X1" "W:X3") are nested in decision list (eg, "W:X1:X10"    "W:X1:X10:X4"), e.g.
                                                            ifelse( N_deci == N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == T, 9,  
                                                                    # SG/INT No. >= truth  & some SG/INT are the same as in the true subgroups;                      & some SG/INT in the decision list (eg "X3<=0 & X1<=0" "X3>0 & X1<=0" "X7<=-0.2 & X1>0" "X7>-0.2 & X1>0") not included in the truth vector (eg "X3<=0 & X1<=0" "X3<=0 & X1>0"  "X3>0 & X1<=0"  "X3>0 & X1>0" )                          
                                                                    ifelse( N_deci == N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) > 0                        & ( length( Deci_vec[Deci_vec %in% Truth_vec]) < length(Deci_vec)  ) , 10, 
                                                                            11 )))))))))
    
  )
  
  #the following is for "subgroup" metric!!!
  sgpct_Nacc_Nalldeci <- length( Deci_vec[Deci_vec %in% Truth_vec])/length(Deci_vec)  #% identified true subgroups among identified subgroups
  sgpct_Nacc_Nalltrue <- length( Deci_vec[Deci_vec %in% Truth_vec])/length(Truth_vec) #% identified true subgroups among true subgroups, 
  
  return( list(value         = value,
               sg_Nacc_Nall  = sgpct_Nacc_Nalldeci , 
               sg_Nacc_NallT = sgpct_Nacc_Nalltrue))
}



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

#e.g. Deci_oneCF_tN1_TF  <- DECI_treeN1_mCF_TF (Deci_oneCF,  vote_1D4b_subgroup, vote_1D3b_subgroup, vote_1D2b_subgroup, vote_1D4b_tree, vote_1D3b_tree, vote_1D2b_tree, tree_true )
#decision=Deci_oneCF
#vote_D4_sg <- vote_1D4b_subgroup
#vote_D3_sg <- vote_1D3b_subgroup
#vote_D2_sg <- vote_1D2b_subgroup
#vote_D4_t  <- vote_1D4b_tree
#vote_D3_t  <- vote_1D3b_tree
#vote_D2_t  <- vote_1D2b_tree
#tree_true <- tree_true

DECI_treeN1_mCF_TF <- function(decision, vote_D4_sg, vote_D3_sg, vote_D2_sg, vote_D4_t, vote_D3_t, vote_D2_t, tree_true){
  #first identify the decision comes from which depth tree 
  if (identical( decision, vote_D4_sg$majority ) ) {
    DECI_treeN1 <<- vote_D4_t$majority[1, "split_variable"] 
  } else if (identical( decision, vote_D3_sg$majority ) ) {
    DECI_treeN1 <<- vote_D3_t$majority[1, "split_variable"] 
  } else if (identical( decision, vote_D2_sg$majority ) ) {
    DECI_treeN1 <<- vote_D2_t$majority[1, "split_variable"] 
  } else if ( identical( decision, "NA" ) ){
    DECI_treeN1 <<- "NA"
  }
  #second check if top node is accurate
  if (identical(tree_true, "NA")==F) {
    TorF= ifelse (identical(DECI_treeN1, tree_true[1, "split_variable"]), 1, 0)
  } else if (tree_true == "NA") {
    TorF= ifelse (identical(DECI_treeN1, tree_true), 1, 0)
  }
  return(TorF)
}

#' Function that check if the tree of the subgroup decision has the correct top node 
#' @param decision the subgroup decision 
#' @param vote_sg 
#' @param vote_t
#' @param tree_true true tree structure
#' 
#' @return subgroup decision   
#' 
#' @export
DECI_treeN1_TREE_TF <- function(decision, vote_sg, vote_t, tree_true){
  if ( identical(decision, "NA")==F ) {
    DECI_treeN1 <<- vote_t[1, "split_variable"] 
  } else if (identical( decision, "NA" ) ){
    DECI_treeN1 <<- "NA"
  }
  
  if ( identical(tree_true, "NA")==F) {
    TorF= ifelse (identical(DECI_treeN1, tree_true[1, "split_variable"]), 1, 0)
  } else if (tree_true == "NA") {
    TorF= ifelse (identical(DECI_treeN1, tree_true), 1, 0)
  }
  return( TorF )
}


#' Function that obtain subgroup decision from other tree-based method such as IT tree   
#' @param vote_subgroup the subgroup decision from voted tree structure
#' 
#' @return subgroup decision   
#' 
#' @export
OTHERTREE_DECI <- function(vote_subgroup){
  if(identical(vote_subgroup,"NA" ) ){
    Deci = "NA" 
  } else {
    Deci <- vote_subgroup$majority
  }
  return(Deci)
}


#' Function that make subgroup defintion as a list of df, each df representing each condition of this subgroup definition and has 3 key variabs: parent_split_var_, parent_sign_, parent_split_val
#' @param Deci_VT_con0 subgroup definition string with "&" sign come from VirtualTWin 
#' 
#' @return subgroup defintion as a list of df, each df representing each condition of this subgroup definition   
#' 
#' @export
SPLITVAR_SIGN_SPLITVAL_LofDf_VT <- function( Deci_VT_con0 ) {
  #step0. remove "&" from original subgroup defintion to prepare
  Deci_VT_con <- stringr::str_split(Deci_VT_con0, " & ", simplify = T) #KEY
  
  #step1. make each subgroup defintion a list, each conditon is an element of the list
  Deci_VT_con_l <- as.list(Deci_VT_con)
  Deci_VT_con_replace1 <- lapply(Deci_VT_con_l,        function(x) stringr::str_replace_all(x, ">=", " >= "))
  Deci_VT_con_replace2 <- lapply(Deci_VT_con_replace1, function(x) stringr::str_replace_all(x, "< ", " < "))
  Deci_VT_con_split    <- lapply(Deci_VT_con_replace2, function(x) stringr::str_split(x, " ")) #KEY
  
  #step2. make a dataframe of this subgroup defintion to manipulate split variable, sign, and split value
  Deci_VT_con_split2 <- data.frame(matrix(unlist(Deci_VT_con_split), nrow=length(Deci_VT_con_split), byrow=TRUE)) %>%
    dplyr::rename(parent_split_var_ = X1, parent_sign_ = X2, parent_split_val_ = X3) %>%
    dplyr::mutate(parent_split_val_ = as.numeric(as.character(parent_split_val_) ) )
  #convert dataframe to list:
  Deci_VT_con_split3 <-  split(Deci_VT_con_split2, seq(nrow(Deci_VT_con_split2)))
  
  return(Deci_VT_con_split3)
}

#' Similar to SPLITVAR_SIGN_SPLITVAL_LofDf_VT, Function that make subgroup defintion as a list of df, each df representing each condition of this subgroup definition and has 3 key variabs: parent_split_var_, parent_sign_, parent_split_val
#' @param Deci_CT_con0 subgroup definition string with "&" sign come from Causal foret or with it's format (VirtualTwin traformed into CF format)
#' 
#' @return subgroup defintion as a list of df, each df representing each condition of this subgroup definition   
#' 
#' @export
SPLITVAR_SIGN_SPLITVAL_LofDf_CT <- function( Deci_CT_con0 ) {
  #step0. remove "&" from original subgroup defintion to prepare
  Deci_CT_con <- stringr::str_split(Deci_CT_con0, " & ", simplify = T) #KEY
  
  #step1. make each subgroup defintion a list, each conditon is an element of the list
  Deci_CT_con_l <- as.list(Deci_CT_con)
  Deci_CT_con_replace1 <- lapply(Deci_CT_con_l,        function(x) stringr::str_replace_all(x, "<=", " <= "))
  Deci_CT_con_replace2 <- lapply(Deci_CT_con_replace1, function(x) stringr::str_replace_all(x, ">", " > "))
  Deci_CT_con_split    <- lapply(Deci_CT_con_replace2, function(x) stringr::str_split(x, " ")) #KEY
  
  #step2. make a dataframe of this subgroup defintion to manipulate split variable, sign, and split value
  Deci_CT_con_split2 <- data.frame(matrix(unlist(Deci_CT_con_split), nrow=length(Deci_CT_con_split), byrow=TRUE)) %>%
    dplyr::rename(parent_split_var_ = X1, parent_sign_ = X2, parent_split_val_ = X3) %>%
    dplyr::mutate(parent_split_val_ = as.numeric(as.character(parent_split_val_) ) )
  #convert dataframe to list:
  Deci_CT_con_split3 <-  split(Deci_CT_con_split2, seq(nrow(Deci_CT_con_split2)))
  
  return(Deci_CT_con_split3)
}


#' Function that transform the splitting format of VirtualTwin style to CausalForest style:
#'1) binary/ordinal splitter >= n.5 to >n, < n.5 to <= n; 2) continuous splitter: >= n to > n-0.001, < n to <= n-0.001
#' @param Deci_VT_con0 subgroup definition string with "&" sign (may come from VirtualTWin or CF)
#' 
#' @return subgroup decision   
#' 
#' @export

VT2CF_format <- function(Deci_VT_con0){
  
  Deci_VT_con_split3 <- SPLITVAR_SIGN_SPLITVAL_LofDf_VT(Deci_VT_con0)
  #step3. adjust splitting format: 1) binary/ordinal splitter >= n.5 to >n, < n.5 to <= n; 2) continuous splitter: >= n to > n-0.001, < n to <= n-0.001
  Deci_VT_con_split4 <- lapply(Deci_VT_con_split3, 
                               function(x) 
                                 if( #scenario 1: binary/ordinal splitter >= 
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) %in% c(2:8) &
                                   x[1,"parent_sign_"]==">=" 
                                 ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = ">", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.5, split_val_round_posi ))
                                 } else if ( 
                                   #scenario 2: binary/oridinal splitter <   
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) %in% c(2:8) &
                                   x[1,"parent_sign_"]     =="<" 
                                 ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = "<=", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.5, split_val_round_posi ))
                                 } else if ( 
                                   #scenario 3: continuous splitter >=
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) > 8 &
                                   x[1,"parent_sign_"]==">=" ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = ">", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.001, split_val_round_posi ))
                                 } else if ( 
                                   #scenario 3: continuous splitter <
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) > 8 &
                                   x[1,"parent_sign_"]=="<" ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = "<=", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.001, split_val_round_posi ))
                                 }
  )
  
  
  #step4. combine information from split variable, sign, and split value into one column 
  Deci_VT_con_split5 <- lapply(Deci_VT_con_split4, function(x) x %>% dplyr::mutate(condition = paste0(parent_split_var_, parent_sign_, parent_split_val_)) %>%
                                 dplyr::select(condition) %>%
                                 `colnames<-`( NULL )
  )
  
  #step5. collapse list of strings into one single string  
  Deci_VT_con_split6 <-  paste ( as.vector( unlist(Deci_VT_con_split5 ) ), collapse = " & " )
  
  return(Deci_VT_con_split6)
  
}


#' Function that convert each subgroup definition into interaction (e.g. a subgroup defined as X1<=0 & X3<=0 will be converted into W:X1:W3)
#' @param condition a specific subgroup definition (one group), e.g.  X1<=0 & X3<=0 
#' @return interaction, e.g. W:X1:W3  
#' 
#' @export

SG_CONDITION_2_INT <- function(condition){
  #make subgroup defintion as a list of df, each df representing each condition of this subgroup definition and has 3 key variabs: parent_split_var_, parent_sign_, parent_split_val
  splitvar_sign_splitval_L <-  SPLITVAR_SIGN_SPLITVAL_LofDf_CT(condition) 
  #combine conditions from a list of a subgroup into a dataframe 
  condition_key_var_df <- dplyr::bind_rows(splitvar_sign_splitval_L, .id = "column_label") 
  #extract key covariates in a vector
  condition_key_var_vec <- unique( str_sort( as.vector( dplyr::pull(condition_key_var_df, parent_split_var_)   ) ) )
  #make it into W:X1:X2 format
  condition_INT <- paste("W", paste(condition_key_var_vec, collapse = ":"), sep = ":")
  
  return(condition_INT)
  
}

#' Function that convert overall subgroup decision into interaction expressions (e.g. subgroup decision as 1: X1<=0; 2: X1>0 will be converted into W:X1 interaction expression)
#' @param condition a specific subgroup definition (one group), e.g.  X1<=0 & X3<=0 
#' @return interaction, e.g. W:X1:W3  
#' 
#' @export
SG2INT <- function(deci){
  if (identical(deci,"NA" )  ) {
    SG2INT_vec = "NA"
  } else if (identical(deci,"Unknown") ) {
    SG2INT_vec = "Unknown"
  } else {
    #extract the subgroup defintion column (e.g. X1<=0, X1>0)
    Deci_col <- deci %>% dplyr::select(subgroup) %>% as.data.frame() %>%  `colnames<-`( NULL ) 
    #make a list containing defintion of each subgroup as each elements of the list
    Deci_col_L <-  split(Deci_col, seq(nrow(Deci_col)))
    #convert each subgroup definition into interaction (e.g. a subgroup defined as X1<=0 & X3<=0 will be converted into W:X1:W3)
    Deci_SG_INT <- lapply(Deci_col_L, SG_CONDITION_2_INT )
    #combine all INT terms into one string, remove duplicate
    SG2INT_str_uniq <- paste(unique(unlist(Deci_SG_INT)), collapse = ' ')
    #turn string into vector
    SG2INT_vec <- strsplit(SG2INT_str_uniq, "\\s+")[[1]]
  }
  return(SG2INT_vec)
}