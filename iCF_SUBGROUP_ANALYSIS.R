#############################################################################################
# subgroup analysis: CATE
# Author: Tiansheng Wang  
# Last update date:12/6/2021
# Version: 0.1         
#############################################################################################
library(tidyverse)
library(broom)

#Personalized package
prop.func <- function(x, trt)
{
  # fit propensity score model
  propens.model <- cv.glmnet(y = trt,
                             x = x, family = "binomial")
  pi.x <- predict(propens.model, s = "lambda.min",
                  newx = x, type = "response")[,1]
  pi.x
}



#-----------------------
#subsetting first!
#-----------------------
#' Function that statify population and saved stratified population (subgroup) in a list. 
#' If decision is from D4, D3, or D2 tree, then it will always split.
#' If decision is from truth, then it could be "NA", i.e. no HTE. In such scenario, set condition_F  = "W<100", i.e. ONE "group" for all observations.  
#' @param decision the final decision obtained from the decision path
#' @param dataset Train or Test dataset
#' 
#' @return the list of stratified population accroding to subgroup decision 
#' 
#' @export

#decision=Deci_D4_oneCFv
#decision=tree_true_subgroup
#dataset=Train
#condition_F= "W <100"

SUBSETTING <- function(decision,  dataset ){
  
  if ( is.null(dim(decision))==T) {
    condition_F= "W <100" 
  } else {
  
  #convert the column of subgroup defintion into a vector of characters for subgroup definition
  condition <- unique(unlist(decision$subgroup))
  #make of list of subgroup defintion
  condition_L <- as.list(condition)
  #replace "& NANANA" condition by empty record to avoid error when subsetting by conditions in the loop step
 # condition_L_noNA <- lapply(condition_L, function(C) ifelse( grepl("& NANANA", C ) , C <- stringr::str_remove(C, " & NANANA"), C) )
  #converto list into a character
  condition_F <- as.character(condition_L)
}

  
  subgrouplist = list()
  for (i in 1:length(condition_F)){
   env=global_env() #Advanced R, P477!!!
    
  #give up eval Fx because EVAL related environment issues can't be fixed (e.g.,  #Error in eval(rlang::parse_expr(condition[i])) : object 'NANANA' not found):
  #solution01: subgroup <- subset(dataset, eval(parse(text = condition[1]) )  )            
  #solution02: subgroup <- subset(dataset, eval(str2lang( condition[1])  )   )                           
  #solution03: subgroup <- subset(dataset, eval( str2expression( condition[i])  )  )  %>% mutate(SubgroupID=i, Definition=condition[i])  
  #solution04: dataset[which(  eval(rlang::parse_expr(condition[i]), env )),]                     
  #other FILTER alternatives:
  #subgroup <- dataset %>% dplyr::filter (  !!str2lang(condition[1])  )  %>% mutate(SubgroupID=i, Definition=condition[i])  
  subgroup <- dataset %>% dplyr::filter (  !!rlang::parse_expr(condition_F[i])  )   %>% dplyr::mutate(SubgroupID=i, Definition=condition_F[i])  
     
  subgrouplist[[i]] <- subgroup
 }
  return(subgrouplist)
  #the 2nd METHOD: this method create a list directly, not need to make an empty list	
# env=global_env() #Advanced R, P477!!!
#  Subgroup <- condition_F %>% purrr::map(~ dataset %>% dplyr::filter_ (.dots = .x) )
#  subKEY_D <- setNames( split(decision, seq(nrow(decision))), rownames(decision))
#  Subgroup_1 = list.zip(subKEY_D, Subgroup)
#  Subgroup_2 <- lapply(Subgroup_1, function(df) cbind(df$subKEY_D, df$Subgroup) )
#  return(Subgroup_2)
}





#' Function that calcualts observed treatment effect size
#' @param subgroup the database (subgroup) in a list of stratified population accroding to subgroup decision
#' 
#' @return The ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) 
#' 
#' @export
#subgroup <- (SUBSETTING (Deci_D4_iCF, Test_ID_m))[[1]]
#SG_D4_y_tr <- lapply(SG_D4_tr, ANA_SUBGROUP)
#subgroup = SG_D4_tr[[1]]

ANA_SUBGROUP <- function(subgroup){
  #----------------------------------------------------------- 
  #won't use forest to predict PS as it seems logistic regression has a better prediction for this simulation setting
  #may consider other method such as SUPER LEARNER later
  
  #regression forest to predict PS
  #---------------------------------------------------------------------------------------------
  #X_s<-subgroup[,c("X1","X2","X3","X4","X5","X6","X8","X9","X10")]
  #Y_s<-subgroup[,"y"]
  #W_s<-subgroup[,"a"]
  #W_s.forest <- grf::regression_forest(X_s, W_s)
  #subgroup$ps <- predict(W_s.forest)$predictions
  
  #-----------------------------------------------------------
  #name <- paste(intTRUE,iterationNo,treeNo, sep="_") #in simulation Fx: Error in paste(intTRUE, iterationNo, treeNo, sep = "_") : 
  env=caller_env()
  #logistic regression to predict PS
  # vars_ex=c("W", "Y", "bl_metformin")
  #W.forest.sg <- grf::regression_forest(subgroup %>% select(-W, -Y), subgroup$W )
  #subgroup$ps <- predict(W.forest.sg)$predictions
 # fmla_ps <- as.formula(paste("W ~ ", paste( colnames(subgroup %>% select (-c(W, Y, SubgroupID, Definition)) ) , collapse= "+")) )
#---------------
#Medicare
#---------------
#change the sign back as "event" in orginal data (in the beginning, I change sign of "event" to make sure the large effect, the better outcome, which followes the pattern of most ML precision methods )
# subgroup$W_f <- factor(subgroup$W, levels = c(0,1), labels=c("SU","DPP4i"))
#subgroup$Y <- ifelse(subgroup$Y==0, 1, 0) 
#subgroup$Y_f <- factor(subgroup$Y, levels=c(0,1), labels=c("MACE", "No MACE"))
  #pt = prop.table(table(Y = subgroup$Y_f, W = subgroup$W_f))
  #_______________________________
  #make a table to get crude RD!
  #reverse order of factor, i.e. + or - sign of RD
  #subgroup %>% dplyr::group_by(W_f, Y_f) %>% dplyr::summarize(n=n()) %>% dplyr::mutate(pct = n / sum(n))
  #subgroup$Y_f = forcats::fct_rev(subgroup$Y_f)

  #fmla_ps <- as.formula(paste("W ~ ", paste(c(vars_forest), collapse= "+")) )
  
# fmla_ps <- as.formula(paste("W ~ ", paste(c( "age_b", "gender"              , "racecat"             , "bl_diabretinopathy"  , "bl_nephropathy" ,      "bl_neuropathy"     ,  "bl_hypertension",     
#                                                        "bl_dyslipidemia" ,     "bl_cerebrovasculardz", "bl_ischemichtdz"     , "bl_peripheralvdz"    , "bl_oedema"      ,      "bl_cardiomyopathy" ,   "bl_arrhythmia" ,      
#                                                        "bl_copd"         ,     "bl_depression"       , "bl_cancer"           , "bl_metformin"        , "bl_glp"         ,      "bl_sglt2"          ,   "bl_lai"       ,       
#                                                        "bl_agi"          ,     "bl_meglitinide"      , "bl_acei"             , "bl_arb"              , "bl_bb"          ,      "bl_ccb"            ,   "bl_statin"    ,       
#                                                        "bl_loop"         ,     "bl_otherdiuretics"   , "bl_hormonet"         , "bl_aspirin"          , "bl_nsaids"      ,      "year"              ,   "uncontroldm"  ,       
#                                                        "hospnodm"        ,     "ervisitdm"           , "phvisit"             , "cvdrvisit"           , "ervisit"        ,      "flushot"           ,   "lics"         ,       
#                                                        "bl_smokingva"     ), collapse= "+")) )
  
 # summary(ps_model <- stats::glm(data = subgroup, fmla_ps , family=binomial("logit")))
 #  subgroup$ps <- predict(ps_model, type = "response")
 # W.RF <- randomForest::randomForest( W ~. , data = subgroup %>% dplyr::select(-Y) )
  #subgroup$ps <- predict(W.forest_sg)$predictions
  #w.out <- WeightIt::weightit(fmla_ps, data = subgroup, estimand = "ATT", method = "ps") # in Simulation Fx: Error: The given response variable, "W", is not a variable in data or the global environment.
  #cobalt::set.cobalt.options(binary = "std")
  #plot occurence of tree
  #png(paste("cov_balance", name, ".png", sep="_")) #In simulation Fx: Error in paste(intTRUE, iterationNo, treeNo, sep = "_") : 
  #cov_balance<-cobalt::love.plot(w.out)
  #print(cov_balance)
  #dev.off()
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ##in Simulation: if leaf size is not appropriate (e.g. too small) then : Error in family$linkfun(mustart) : Argument mu must be a nonempty numeric vector
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #Probability of receving Tx, as a function of covatiates
  #summary(ps_model <- glm(data = subgroup, fmla_ps, family=binomial("logit"))) 
  #Use the model to predict the probability of Tx for each observation
  #subgroup$ps <- predict(ps_model, type = "response")
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
  #Plot propensity score 
  #png(paste("PS_distrib", name, ".png", sep="_")) #in simulation Fx: Error in paste(intTRUE, iterationNo, treeNo, sep = "_") : 
  #PS.distribution<-ggplot(data=subgroup,aes(x=ps, group=factor(W_f), fill=factor(W_f)))+
  #  geom_histogram(aes(y=..density..),alpha = 0.75,binwidth=0.004,position = position_dodge(width=0.01))+
  #  theme_classic()+
  #  xlab(paste("Predicted probability of Receiving Treatment", sep=""))+ylab("Density")+
  #  labs(fill = "Observed")
  #print(PS.distribution)
  #dev.off()
  #Alternative: use IPTW to estimate weighted risks and measures of effect
  subgroup$iptw_u <- ifelse(subgroup$W==1, 1/subgroup$ps, 1/(1-subgroup$ps)) 
#  subgroup$iptw_tsf <- ifelse(subgroup$W==1, 1/subgroup$ps, -1/(1-subgroup$ps)) 
#  subgroup$Y_star <-  subgroup$Y*iptw_tsf
  
    p_exposure <- sum(subgroup$W) / nrow(subgroup)
  subgroup$iptw_s <- ifelse(subgroup$W==1, p_exposure/subgroup$ps, (1-p_exposure)/(1-subgroup$ps))
  subgroup$smrw <- ifelse(subgroup$W==1, 1, (subgroup$ps)/(1-subgroup$ps))
  ####################################################################################

  DLT_Y <- function(subdataset, Wt){
    
  if (length(unique(Y)) >=8){
      
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
  } else if (length(unique(Y)) ==2){
    # starting value helps: https://stackoverflow.com/questions/35618026/what-do-these-r-glm-error-messages-mean-error-no-valid-set-of-coefficients-ha
    # length of 'start' should equal 2 and correspond to initial coefs for c("(Intercept)", "W")
    if(Wt == "iptw_s"){ dltY_Mod<- spaMM::spaMM_glm(Y ~ W,
                                             data=subdataset, 
                                             family= binomial(link="identity"),  #risk difference scale, 
                                             weights=iptw_s#, 
                                              )
    } else if (Wt == "smrw") {  dltY_Mod <- spaMM::spaMM_glm(Y ~  W,
                                                      data=subdataset, 
                                                      family= binomial(link="identity"), #risk difference scale
                                                      weights=smrw#,
                                                      ) 
    } else if (Wt == "none") { dltY_Mod <- spaMM::spaMM_glm(Y ~  W,
                                                     data=subdataset, 
                                                     family= binomial(link="identity")#, #risk difference scale
                                                    ) 
    }
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
#TorChi <- t.test(x=subgroup$Y[subgroup$W==1], y = subgroup$Y[subgroup$W==0],
#       alternative = "two.sided",
#       mu = 0, paired = FALSE, var.equal = FALSE,
#       conf.level = 0.95
#       )
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



#' Function that extract calcualts observed treatment effect size from each depth of iCF and saved in a separate dataframe
#' In order to record data from each simulation
#' @param Subgroup_dltY_D the list of stratified population accroding to subgroup decision
#' 
#' @return The dataframe of ATE (IPTW), ATT (SMR), and crude treatment effect (i.e., Delta Y) 
#' 
#' @export
#' 

#EXTRACT_SUB_ANA(sg_D4_y_tr_un,4,"train", methodL)
#Subgroup_dltY_D=sg_D4_y_tr_un  
#Depth=4 
#datatype="train" 
#method=methodL

EXTRACT_SUB_ANA <- function(Subgroup_dltY_D,  Depth, datatype, method) {
  if (datatype=="train") {
                          data_label="tr"
                          } else if (datatype=="test") {
                            data_label="te"
                          }
  #Extract info from dataframes in the Subgroup_dltY list 
  SubID         <- unlist(lapply(Subgroup_dltY_D,function(x) x$SubgroupID[1] ) )
  SubDef        <- unlist(lapply(Subgroup_dltY_D,function(x) x$Definition[1])) 
  
  
#  SubID         <- unlist(lapply(Subgroup_dltY_D,function(x) x[1,13]))
#  SubDef        <- unlist(lapply(Subgroup_dltY_D,function(x) x[1,14])) 

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
#' 
#' 
#' 

#Deci_D4=vote_D4_subgroup.ori$majority
#Deci_D3=vote_D3_subgroup.ori$majority
#Deci_D2=vote_D2_subgroup.ori$majority
#traindata=Train
#testdata=Test_ID
#methodL="iCF"

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
#' @param Deci_IT the decision from IT 
#' @train Training dataset
#' @test Testing dataset
#' @methodL method label, "IT"
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



#' Function that obtain residual and MSE for each observation 
#' @param decision the subgroup decision by different methods 
#' @param dataset the dataset analyzed for residual & MSE
#' @param label if residual and MSE is for true effect or observed effect
#' 
#' @return the dataframe of residual & MSE
#' 
#' @export
#' 

DltY_DATA <- function (decision, dataset, label){
  #generate a list of subgroups (i.e.stratified population by subgroup definitions)
  subgroup_L        <- SUBSETTING(decision, dataset)
  #obtain treatment effect (ATE, ATT) of each subgroup
  DltY_L            <- lapply(subgroup_L, ANA_SUBGROUP)
  #remove useless iterms of each sublist
  DltY_L_l_dataonly <- lapply(DltY_L, function(x) rlist::list.remove(x, c("model_ate", "model_att", "model_cru", "TorChi")) )
  #combine all obseravtions after removing useless iterms 
  DltY_L_dataonly   <- lapply(DltY_L_l_dataonly, function(x) data.table::rbindlist (x)  )
  
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


#' Function that obtain subgroup ID 
#' @param decision the subgroup decision by different methods 
#' @param dataset the dataset analyzed for residual & MSE
#' 
#' @return the dataframe of residual & MSE
#' 
#' @export
#' 

GET_SUBGROUP_ID <- function (decision, dataset){
  #generate a list of subgroups (i.e.stratified population by subgroup definitions)
  subgroup_L        <- SUBSETTING(decision, dataset)

  DltY <- rlist::list.rbind( subgroup_L ) %>% #combine all observations from each subgroups into one dataset 
    dplyr::arrange (ID)                    #order by ID
  
  return(DltY)
}

#' Function that extract MSE from MSE dataframe statify population and saved stratified population (subgroup) in a list. 
#' @param DltY the dataframe of MSE obtained from function DltY_DATA
#' 
#' @return the list of MSE, residual by IPTW (ate) or SMR (att) weighting
#' 
#' @export
#' 

MSE_DltY <- function(DltY_true, DltY) {
  env=global_env()
  MSE_df <- dplyr::left_join(DltY_true, dplyr::select(DltY,  c(ID, SubgroupID, Definition,     ps,            iptw_u,     iptw_s,        smrw,
                                                               dltY_ate, dltY_ate_low, dltY_ate_up, dltY_att, dltY_att_low, dltY_att_up, dltY_crude, dltY_crude_low, dltY_crude_up)), 
                             by = "ID" ) %>%
            dplyr::filter(SubgroupID.y != "NA") %>% #when left_join with decision by aVirtualTwin, the total sample size < N (some population not included in the final decision, which leads to )
    dplyr::mutate(residual_ate =  dltY_ate_t -  dltY_ate,
                   mse_ate = mean(residual_ate^2),
                   residual_att =  dltY_att_t -  dltY_att,
                   mse_att = mean(residual_att^2))
  
  return(list(residual_ate  = MSE_df[1, "residual_ate"] , 
              MSE_ate       = MSE_df[1, "mse_ate"], 
              residual_att  = MSE_df[1, "residual_att"] , 
              MSE_att       = MSE_df[1, "mse_att"])
  )
}

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

#' Function that identify subgroup with the largest CATE defined by: 
#' 1) binary outcome: Q(A)=(P(Y=1|T=1, X∈A) - P(Y=1|T=0, X∈A)) – (P(Y=1|T=1) - P(Y=1|T=0))
#' 2) continuous:     Q(A)=(E(Y=1|T=1, X∈A) - E(Y=1|T=0, X∈A)) – (E(Y=1|T=1) - E(Y=1|T=0))
#' @param DltY.data calculated dltY
#' 
#' @return the subgroup leading to the largest Q(A)
#' @export
#' 
CATE_MAX_SG <- function(DltY.data){
  
  CATE_L_true_te  <-   DltY.data %>%
    dplyr::select(SubgroupID, Definition, dltY_ate) %>%
    group_by(Definition) %>%
    summarise(dltY_ate = mean(dltY_ate), n=n()) %>%
    dplyr::mutate(Dlt_ATEvsCATE = (dltY_ate - ATE_all_te) ) %>%
    dplyr::mutate(Dlt_ATEvsCATE_abs = abs(dltY_ate - ATE_all_te) ) %>%
    dplyr::arrange(Dlt_ATEvsCATE_abs)
  
  CATE_largestSG_true_te <- CATE_L_true_te$Definition[1]
  
  if ( CATE_largestSG_true_te=="W <100") {CATE_largestSG_true_te = "NA"}
  
  return(CATE_largestSG_true_te)
  
}
