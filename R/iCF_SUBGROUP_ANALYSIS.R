#' prop.func
#' 
#' Function that predict PS using LASSO algorithm, from the Personalized package
#' @param x covariates
#' @param trt "treatment" 
#'  
#' @return the prediced propensity score
#' 
#' @export

prop.func <- function(x, trt)
{
  # fit propensity score model
  propens.model <- cv.glmnet(y = trt,
                             x = x, family = "binomial")
  pi.x <- predict(propens.model, s = "lambda.min",
                  newx = x, type = "response")[,1]
  pi.x
}


#' GET_SUBGROUP_ID
#' 
#' Function that obtain subgroup ID 
#' @param decision the subgroup decision by different methods 
#' @param dataset the dataset analyzed for residual & MSE
#' 
#' @return the dataframe of residual & MSE
#' 
#' @export

GET_SUBGROUP_ID <- function (decision, dataset){
  #generate a list of subgroups (i.e.stratified population by subgroup definitions)
  subgroup_L        <- SUBSETTING(decision, dataset)
  
  DltY <- rlist::list.rbind( subgroup_L ) %>% #combine all observations from each subgroups into one dataset 
    dplyr::arrange (ID)                    #order by ID
  
  return(DltY)
}

#-----------------------
#subsetting first!
#-----------------------
#' SUBSETTING
#' 
#' Function that statify population and saved stratified population (subgroup) in a list. 
#' If decision is from D5, D4, D3, or D2 tree, then it will always split.
#' If decision is from truth, then it could be "NA", i.e. no HTE. In such scenario, set condition_F  = "W<100", i.e. ONE "group" for all observations.  
#' @param decision the final decision obtained from the decision path
#' @param dataset Train or Test dataset
#' 
#' @return the list of stratified population accroding to subgroup decision 
#' 
#' @export

SUBSETTING <- function(decision,  dataset ){
  if ( is.null(dim(decision))==T) {
    condition_F= "W <100" 
  } else {
  #convert the column of subgroup defintion into a vector of characters for subgroup definition
  condition <- unique(unlist(decision$subgroup))
  #make of list of subgroup defintion
  condition_L <- as.list(condition)
  condition_F <- as.character(condition_L)
}
  subgrouplist = list()
  for (i in 1:length(condition_F)){
   env=global_env() #Advanced R, P477!!!
  subgroup <- dataset %>% dplyr::filter (  !!rlang::parse_expr(condition_F[i])  )   %>% dplyr::mutate(SubgroupID=i, Definition=condition_F[i])  
     
  subgrouplist[[i]] <- subgroup
 }
  return(subgrouplist)
}






#' ANA_SUBGROUP
#' 
#' Function that calcualts observed treatment effect size
#' @param subgroup the database (subgroup) in a list of stratified population accroding to subgroup decision
#' 
#' @return The ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) 
#' 
#' @export

ANA_SUBGROUP <- function(subgroup){
  #----------------------------------------------------------- 
  #may consider other method such as regression_foest, SUPER LEARNER later
  env=caller_env()
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ##in Simulation: if leaf size is not appropriate (e.g. too small) then : Error in family$linkfun(mustart) : Argument must be a nonempty numeric vector
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------

  #------------------------------------------
  # LOGISTIC REGRESSION w/ LASSOP PENALTY
  #------------------------------------------
  if ("ID" %in% colnames(subgroup)) {
  subgroup$ps <- prop.func(as.matrix( subgroup %>% dplyr::select (-c(W, Y, SubgroupID, Definition, ID))), #matrix of covariates
                           subgroup$W #outcome for glmnet
                           ) 
  } else {
  subgroup$ps <- prop.func(as.matrix( subgroup %>% dplyr::select (-c(W, Y, SubgroupID, Definition))), #matrix of covariates
                             subgroup$W #outcome for glmnet
    )  
    
  }
  #Distribution of propensity scoreE87AA3DA|...
  summary(subgroup$ps[subgroup$W==1]) #...among exposued 
  summary(subgroup$ps[subgroup$W==0]) #...among unexposed 

  subgroup$iptw_u <- ifelse(subgroup$W==1, 1/subgroup$ps, 1/(1-subgroup$ps)) 
  
    p_exposure <- sum(subgroup$W) / nrow(subgroup)
  subgroup$iptw_s <- ifelse(subgroup$W==1, p_exposure/subgroup$ps, (1-p_exposure)/(1-subgroup$ps))
  subgroup$smrw <- ifelse(subgroup$W==1, 1, (subgroup$ps)/(1-subgroup$ps))
  ####################################################################################

  #' DLT_Y
  #' 
  #' Function that calcualts observed treatment effect size, leave this function HERE!!!
  #' @param subdataset the database (subgroup) in a list of stratified population accroding to subgroup decision
  #' @param Wt weight 
  #' 
  #' @return The ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) 
  #' 
  #' @export
  
  DLT_Y <- function(subdataset, Wt){
    
   # if (length(unique(Y)) >=8){
      
      if (Wt == "iptw_s"){dltY_Mod <- glm(Y ~ W,
                                          data=subdataset, 
                                          family=  gaussian,  
                                          weights=iptw_s)
      } else if (Wt == "smrw") {  dltY_Mod <- glm(Y ~  W,
                                                  data=subdataset, 
                                                  family= gaussian, 
                                                  weights=smrw) 
      } else if (Wt == "none") {  dltY_Mod<- glm(Y ~  W,
                                                 data=subdataset, 
                                                 family= gaussian) 
      } 
    
    #exposure coefficient, i,e.the difference in the expected value of Y between the exposed and unexposed
    dlt_Y <- summary(dltY_Mod)$coefficient
    dltY_CI  <- confint(dltY_Mod)
    dltY     <-  cbind.data.frame(dlt_Y,  dltY_CI) [2,] %>% dplyr::mutate(Sub_def= subgroup$Definition[1] ) %>% dplyr::rename(conf.low="2.5 %", conf.high="97.5 %")
    
    return(dltY)
    
  }
  
    
  dltY_ate <- DLT_Y(subgroup, "iptw_s")
  dltY_att <- DLT_Y(subgroup, "smrw")
  dltY_crude <- DLT_Y(subgroup, "none")
#dplyr::mutate doesn't work, can't add these in oneline
  #ate
subgroup$dltY_ate      <- round(  dltY_ate$Estimate  , digits=2)
subgroup$dltY_ate_low  <- round(  dltY_ate$conf.low  , digits=2)
subgroup$dltY_ate_up   <- round(  dltY_ate$conf.high  , digits=2)
  #att
subgroup$dltY_att      <- round(  dltY_att$Estimate  , digits=2)
subgroup$dltY_att_low  <- round(  dltY_att$conf.low  , digits=2)
subgroup$dltY_att_up   <- round(  dltY_att$conf.high  , digits=2)
  #crude
subgroup$dltY_crude     <- round(  dltY_crude$Estimate , digits=2)
subgroup$dltY_crude_low <- round(  dltY_crude$conf.low  , digits=2)
subgroup$dltY_crude_up  <- round(  dltY_crude$conf.high  , digits=2)

if (length(unique(Y)) >=8){

TorChi <- t.test(subgroup$Y~ subgroup$W)
subgroup$TorChi = round(TorChi$statistic,3)
subgroup$TorChi_df = TorChi$parameter
subgroup$TorChi_pval = round(TorChi$p.value,5)

} else if  (length(unique(Y)) ==2){
TorChi <- chisq.test( table(subgroup$Y, subgroup$W), correct = FALSE)
 subgroup$TorChi = round(TorChi$statistic,3)
 subgroup$TorChi_df = TorChi$parameter
 subgroup$TorChi_pval = round(TorChi$p.value,5)
}

  return(list(subgroup=subgroup,
              model_ate=dltY_ate,
              model_att=dltY_att,
              model_cru=dltY_crude,
              TorChi=TorChi))
}



#' EXTRACT_SUB_ANA
#' 
#' Function that extract calculates observed treatment effect size from each depth of iCF and saved in a separate dataframe
#' In order to record data from each simulation
#' @param Subgroup_dltY_D the list of stratified population accroding to subgroup decision
#' 
#' @return The dataframe of ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) 
#' 
#' @export

EXTRACT_SUB_ANA <- function(Subgroup_dltY_D,  Depth, datatype, method) {
  if (datatype=="train") {
                          data_label="tr"
                          } else if (datatype=="test") {
                            data_label="te"
                          }
  #Extract info from dataframes in the Subgroup_dltY list 
  SubID         <- unlist(lapply(Subgroup_dltY_D,function(x) x$SubgroupID[1] ) )
  SubDef        <- unlist(lapply(Subgroup_dltY_D,function(x) x$Definition[1])) 

assign ( paste("dltY_ate",      data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_ate[1])) )
assign ( paste("dltY_ate_low",  data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_ate_low[1])) )
assign ( paste("dltY_ate_up",   data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_ate_up[1])) )
assign ( paste("dltY_att",      data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_att[1])) )
assign ( paste("dltY_att_low",  data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_att_low[1])) )
assign ( paste("dltY_att_up",   data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_att_up[1])) )
assign ( paste("dltY_crude",    data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_crude[1])) )
assign ( paste("dltY_crude_low",data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_crude_low[1])) )
assign ( paste("dltY_crude_up", data_label,sep = "_")  ,  unlist(lapply(Subgroup_dltY_D,function(x) x$dltY_crude_up[1])) )

assign ( paste("TorChi",    data_label,sep = "_") ,  unlist(lapply(Subgroup_dltY_D,function(x) x$TorChi[1])) )
assign ( paste("TorChi_df",   data_label,sep = "_") ,  unlist(lapply(Subgroup_dltY_D,function(x) x$TorChi_df[1])) )
assign ( paste("TorChi_pval", data_label,sep = "_") ,  unlist(lapply(Subgroup_dltY_D,function(x) x$TorChi_pval[1])) )

#set environment, #Advanced R, P477!!!
env=global_env() 
#make a dataframe: use cbind.data.frame to avoid rbind()/cbind() conversion from numeric to factor!!!
DltY_D <- as.data.frame(cbind.data.frame(SubID, SubDef, eval(rlang::parse_expr(paste("dltY_ate",      data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_ate_low",   data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_ate_up",    data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_att",       data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_att_low",   data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_att_up",    data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_crude",     data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_crude_low", data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("dltY_crude_up",  data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("TorChi",         data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("TorChi_df",      data_label,sep = "_") ) ),
                                                       eval(rlang::parse_expr(paste("TorChi_pval",    data_label,sep = "_") ) )
                                             ) 
                                   )
#assign name to the above generated dataframe
names(DltY_D) <- c(paste("SubID",                     method, sep = "_")  , 
                   paste("SubDef",                    method, sep = "_")  ,
                   paste("dltY_ate",      data_label, method, sep = "_")  , 
                   paste("dltY_ate_low",  data_label, method, sep = "_")  , 
                   paste("dltY_ate_up",   data_label, method, sep = "_")  , 
                   paste("dltY_att",      data_label, method, sep = "_")  ,  
                   paste("dltY_att_low",  data_label, method, sep = "_")  ,  
                   paste("dltY_att_up",   data_label, method, sep = "_")  ,  
                   paste("dltY_crude",    data_label, method, sep = "_")  ,  
                   paste("dltY_crude_low",data_label, method, sep = "_")  , 
                   paste("dltY_crude_up", data_label, method, sep = "_")  ,
                   paste("TorChi",        data_label, method, sep = "_")  ,
                   paste("TorChi_df",     data_label, method, sep = "_")  , 
                   paste("TorChi_pval",   data_label, method, sep = "_")  
                   ) 
#add depth info
DltY_D[[paste("Depth", method, sep = "_")]] <- Depth
return(DltY_D)
}




#' DltY_DATA
#' 
#' Function that obtain residual and MSE for each observation, used in CATE_SG function below 
#' @param decision the subgroup decision by different methods 
#' @param dataset the dataset analyzed for residual & MSE
#' @param label if residual and MSE is for true effect or observed effect
#' 
#' @return the dataframe of residual & MSE
#' 
#' @export

DltY_DATA <- function (decision, dataset, label){
  #generate a list of subgroups (i.e.stratified population by subgroup definitions)
  subgroup_L        <- SUBSETTING(decision, dataset)
  #obtain treatment effect (ATE, ATT) of each subgroup
  DltY_L            <- lapply(subgroup_L, ANA_SUBGROUP)
  #remove useless iterms of each sublist
  DltY_L_l_dataonly <- lapply(DltY_L, function(x) rlist::list.remove(x, c("model_ate", "model_att", "model_cru", "TorChi")) )
  #combine all obseravtions after removing useless iterms 
  
 #remvoe data.table function  by Tian 9/19/2023 
 #DltY_L_dataonly   <- lapply(DltY_L_l_dataonly, function(x) data.table::rbindlist (x)  )
  
  
  DltY_L_dataonly   <- lapply(DltY_L_l_dataonly, function(x) dplyr::bind_rows(x) )
  

  oldnames = c("dltY_ate",   "dltY_ate_low",   "dltY_ate_up",   "dltY_att",   "dltY_att_low",   "dltY_att_up",   "dltY_crude",   "dltY_crude_low",   "dltY_crude_up", "TorChi", "TorChi_df", "TorChi_pval")
  
  newnames = c("dltY_ate_t", "dltY_ate_low_t", "dltY_ate_up_t", "dltY_att_t", "dltY_att_low_t", "dltY_att_up_t", "dltY_crude_t", "dltY_crude_low_t", "dltY_crude_up_t", "TorChi_t", "TorChi_df_t", "TorChi_pval_t")
  
  
  DltY <- rlist::list.rbind(DltY_L_dataonly) %>% #combine all observations from each subgroups into one dataset 
          dplyr::arrange (ID)                    #order by ID
    
  if (label == "truth"){
    DltY2 <- DltY %>%                      #order by ID
      dplyr::rename_at(vars( all_of(oldnames)), ~ newnames) #rename variable for dataframe of truth to prepare for data linkage
  } else {
    DltY2 <- DltY
  }
  return(DltY2)
}


#' CATE_SG
#' 
#' Function that calculate per-subgroup treatment effect: 
#' @param SGdecision subgroup decision
#' @param dat data set for overall population
#' 
#' @return the subgroup leading to the largest Q(A)
#' @export

CATE_SG <- function(dat, SGdecision){
  ID <-1:nrow(dat)
  dat_ID <- cbind(dat, as.vector(ID)) %>% dplyr::rename (ID=`as.vector(ID)`) 
  Delta_Y <- DltY_DATA(SGdecision, dat_ID, "non-test")
  
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


