####################################################################
# SUBGROUP DECISION 
# Author: Tiansheng Wang
# Last update date: 12/8/2020
# Version: 0.1      
####################################################################

#' TREE2SUBGROUP
#' 
#' Wrapper function that convert tree structure to subgroup 
#' @param tree tree structure df
#'  
#' @return subgroup in list format
#' 
#' @export

TREE2SUBGROUP <- function(tree){
  subgroup.ori    <- PRE_MAJORITY_SUBGROUP ( list(tree ) )
  
  subgroup_key.ori   <-lapply(subgroup.ori,   function(df) subset(df, select=c("subgroupID", "subgroup")))
  
  #did not vote for subgroup, use the synthetic voted tree structure (i.e. averaging splitting value) to get subgroup decision
  vote_subgroup.ori <- list(majority=subgroup_key.ori[[1]])
  
  return(vote_subgroup.ori )
}



#' SGMODEL_DATA
#' 
#' Function that does prepare data that:
#' 1) ready to apply the formulas including interacting terms between treatment and factor subgroup deicision variable; 
#' 2) transformed outcome
#' @param dat data set
#' @param outcome_type "actual" or "transform" 
#'https://gsbdbi.github.io/ml_tutorial/hte_tutorial/hte_tutorial.html
#' Estimation of Heterogeneous Treatment Effects - prepared for “Machine Learning and Causal Inference” class
#'  
#' @return the prepared data and "contr"
#' 
#' @export


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

#' PICK_m2345 
#' 
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


#' PICK_CV
#' 
#' Function that does shows the depth of subgroup selected 
#' @param HTE_P_value P-value for raw full CF
#' @param MAF_m   the mean MSE value across CV for main effect model 
#' @param MAF_g2  the mean MSE value across CV for D2 SG model   
#' @param MAF_g3  the mean MSE value across CV for D3 SG model 
#' @param MAF_g4  the mean MSE value across CV for D4 SG model 
#' @param MAF_g5  the mean MSE value across CV for D5 SG model
#' @param V_D2_SG.L  the list of voted subgroup at D2  
#' @param V_D3_SG.L  the list of voted subgroup at D3 
#' @param V_D4_SG.L  the list of voted subgroup at D4 
#' @param V_D5_SG.L  the list of voted subgroup at D5
#' @param S_D2_T_r  the stability of relaxed (ignoring splitting values) of voted tree structure at D2   
#' @param S_D3_T_r  the stability of relaxed (ignoring splitting values) of voted tree structure at D3 
#' @param S_D4_T_r  the stability of relaxed (ignoring splitting values) of voted tree structure at D4 
#' @param S_D5_T_r  the stability of relaxed (ignoring splitting values) of voted tree structure at D5
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









#' CF_GROUP_DECISION
#' 
#' Function that obtain models from subgroup decisions at different depth
#' @param HE_P_value raw heterogeneity test P-value 
#' @param dat the dataset 
#' @param method_CF method for CF: "iCF" or "oneCF"
#' @param outcome_type_at "actual" or "transform"
#' @param P_threshold threshold for P-value of heterogeneity  
#'  
#' @return subgroup decision
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



