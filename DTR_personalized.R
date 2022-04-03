

EXTRACT_INT <- function(SI_results_df){
  
if (nrow(SI_results_df) >=1) {
  INT <- dplyr::pull(SI_results_df, variable_int) 
} else {
  INT <- "NA"
}
return(INT)
}

#DynTxRegime package for RWL & OWL


## 

#' Function that define value function to find the optimal treatment plan for individuals undr the single-decision time setting
#' @param y the actual outcome
#' @param assignment the treatment assignment
#' @param decision the optimal treatment
#' @param propensity the group level marginal propensity score
#' 
#' @return the list of stratified population accroding to subgroup decision 
#' 
#' @export
value = function(y, assignment, decision, propensity) {
  v = sum(y * (assignment == decision) / propensity)
  match.fac = sum((assignment == decision) / propensity)
  return(v/match.fac)
}


DynTxRegime_WL <- function(method_WL, method_kernel, lambdas, cvFolds){
  
  WL.decision=rep(0,nrow(Train) )
  WL.value=c()
  K=10
  foldid = sample(1:K, size = nrow(Train), replace = TRUE)
  
  for (k in 1:K){
    idx = which(foldid==k)
    train = subset(Train, foldid != k)
    test  = subset(Train, foldid == k)
    if (method_WL=="RWL"){
      fit.WL <- DynTxRegime::rwl(moPropen = moPropen, 
                                  moMain = moMain,
                                  ## define a linear model for the regime funtion
                                  regime = as.formula(paste0(" ~ ", paste0(colnames(X), collapse = " + "))),
                                  data = train,
                                  reward = train$Y, txName = "W",
                                  kernel = method_kernel, ## choose the linear kernel
                                  lambdas = lambdas, ## penalty tuning parameter
                                  cvFolds = cvFolds, ## for selecting the tuning parameters
                                  verbose = T #suppress screen prints
                                  )
    } else if (method_WL=="OWL"){
      fit.WL <- DynTxRegime::owl(moPropen = moPropen, 
                                  regime =  as.formula(paste0(" ~ ", paste0(colnames(X), collapse = " + "))),
                                  data = train, 
                                  reward = train$Y,  
                                  txName = "W",
                                  kernel = method_kernel,
                                  lambdas = lambdas, ## penalty tuning parameter
                                  cvFolds = cvFolds,
                                  verbose = F
                                 )  
    } else if (method_WL =="AIPWE"){
      fit.WL <- DynTxRegime::optimalSeq(moPropen = moPropen, 
                                        moMain = moMain, 
                                        moCont = moCont, 
                                        regimes = regime_aipwe,
                                         data = train, 
                                         response = train$Y, 
                                         txName = "W",
                                         Domains = cbind(rep(-1,ncol(X)+1),rep(1,ncol(X)+1)), ## domain for eta
                                         pop.size = nrow(train), 
                                         starting.values = rnorm(ncol(X)+1, 0, 0.1)
                                        )
    }
    #find the optimal treatment for test data
    WL.decision[idx] = DynTxRegime::optTx(x=fit.WL, newdata = test)$optimalTx
    WL.value = c(WL.value, value(test$Y, Train[idx, "W"], WL.decision[idx], prop[idx]))
  }
  
  Train_opt <- cbind(Train, WL.decision) %>% rename(optTx = WL.decision)
  
  
  if ( length(unique( Train_opt$optTx) ) == 1 #| min(as.vector(table(Train_opt$optTx))) < nrow(Train_opt)/20, won't require sampel size for now
       ) {
    cov_int_df="NA"
    INT ="NA"
  } else if ( length(unique( Train_opt$optTx) ) > 1 ) {
    
    #will figure out automatically identify continuous and discrete variables later.
    #Train_opt_cont <- Train_opt %>% select(X2, X4, X7, X10, optTx)
    #Train_opt_disc <- Train_opt %>% select(X1, X3, X5, X6, X8, X9, optTx)
    long_opt_cont <- Train_opt %>% 
                      dplyr::select(X2, X4, X7, X10, optTx) %>% 
                      tidyr::pivot_longer(-optTx, names_to = "variables", values_to = "value")
    long_opt_disc <- Train_opt %>% 
                     dplyr::select(X1, X3, X5, X6, X8, X9, optTx) %>%
                     tidyr::pivot_longer(-optTx, names_to = "variables", values_to = "value")
    #_________________________________________
    #t-test for continuous variables
    stat.test_cont <- long_opt_cont %>%
                      group_by(variables) %>%
                      rstatix::t_test(value ~ optTx) 
    #_________________________________________
    #chi-square test for discrete variables
    #group by variables
    vnames_disc = c("X1", "X3", "X5","X6","X8", "X9")
    vnames_cont = setdiff(colnames(X) , vnames_disc)
    #long_opt_cont <- Train_opt_disc %>% select(-optTx)
    var_list <- list()
    for (v in 1 : length(vnames_cont)){
      var <- long_opt_disc %>% filter (variables== vnames_disc[v])
      var_list[[v]]=var
    }
    #chi-square test
    var_list_chi <- lapply(var_list, function(x)  rstatix::chisq_test(as.table( table (x$value, x$optTx ) ) ) %>% mutate(variables = x$variables[1] ) )
    stat.test_disc <- do.call("rbind", var_list_chi)%>% select(variables, everything())
    
    #combine t-test and chi-square test p-value for Hommel adjustment
    P_Hommel <-rbind(stat.test_cont %>% select(-n1, -n2, -group1, -group2, -`.y.`) %>% mutate(method="t-test"), 
                     stat.test_disc %>% select (-n, -p.signif))  %>%
      rstatix::adjust_pvalue(method = "hommel") %>%
      rstatix::add_significance()
    cov_int_df <- P_Hommel %>%
      filter(p.adj.signif != "ns", 
             p.adj <  Hommel_P_threshold) %>%
      dplyr::arrange ( p.adj) %>%
      dplyr::mutate (variable_int = paste0("W:", variables))
    
    INT <- EXTRACT_INT(cov_int_df )
    
  }
  
  return(list(INT=INT,
              cov_int_df = cov_int_df,
              WL.decision = WL.decision,
              WL.value = WL.value))
}




# ----------------------------------------------------
# AIPWE  - Baqun Zhang ..., very time-consuming, give up!
# ----------------------------------------------------
# Regime function as defined in the paper
#regime_aipwe <<- function(eta0, eta1, eta2, eta3, eta4, eta5, eta6,
#                   eta7, eta8, eta9, eta10, data) {
#  eta = matrix(c(eta0, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10))
#  X = cbind(1, data$X1, data$X2, data$X3, data$X4, data$X5, data$X6, data$X7, data$X8, data$X9, data$X10)
#  tst <- c(X %*% eta > 0)
#  rec <- rep(1, nrow(x = data))
#  rec[!tst] <- 0
#  return( rec )
#}
#s_aipwe =  Sys.time()
#AIPWE <- DynTxRegime_WL(method_WL = "AIPWE")
#INT_OWL=OWL$INT
#e_aipwe =  Sys.time()
#time_aipwe = e_aipwe - s_aipwe

#Deci_stima <- OTHERTREE_DECI(vote_stimaDynTxRegime::propen(object = result)_subgroup)






#PSOPTMATCH <- function(dat0){
#  psform1 <- as.formula(paste0("a ~ ", paste0("X", 1:10, collapse = " + ")))
#  aGlm <- glm(psform1, family=binomial(), data=dat0)
#  ps <- predict(object = aGlm, type = "response")
#  ps_m<-match_on(aGlm, caliper = 0.25*sd(ps))
#  fm <- fullmatch(ps_m, data = dat0)
#  dat_m <- na.omit(cbind(dat0, fm))
#  return(dat_m)
#}



#wont' plot as these are by qualitative subgroup (recommended vs not recommended)
#plot(subgrp.model)
#plot(subgrp.model, type="density")
#plot(subgrp.model, type="interaction")
#plot(subgrp.model, type="conditional")
#=============================
#sg_mod_row <- fit.subgroup(x = as.matrix(X_a), y = Y, trt = W_a,
#                             propensity.func = prop.func,
#                             method="weighting",
#                             loss = "owl_logistic_loss_lasso",
#                             nfolds = 10)
#comp_row <- summarize.subgroups(sg_mod_row)
#cov_int_row_df<- as.data.frame(print(comp_row, p.value = 1e-9))
#cov_int_row <- rownames(cov_int_row_df)




#preapre dataset format for SIDES method
#Train_m_sides <- Train_m %>% select(-c("z.a_trueps", "Y1", "fm")) %>%
#                 select(a, everything()) %>%
#                 select(Y, everything())

#s1 = SIDES(all_set=Train_m %>% select(-c("z.a_trueps", "Y1", "fm")) %>%
#             select(a, everything()) %>%
#             select(Y, everything()) #the "global data set", the 1st column =outcome, the 2nd column = Tx, others=covariates
#          ,
#           type_var=c(rep("continuous",4),rep("ordinal",6)), 
#           type_outcome="continuous",
#           level_control=0, D=0, L=3, S=100, M=5, gamma=c(1,1,1), H=1, num_crit=1,
#           alpha=0.10, nsim=1000, ord.bin=10, upper_best=TRUE, seed=42)


#assessing stima: 

#if(vote_stima_subgroup=="NA"){
#  Deci_stima= "NA" 
#} else {
#  Deci_stima <- vote_stima_subgroup$majority
#}
#check if true

#TF_D_stima <- ifelse(identical(Deci_stima,  tree_true_subgroup), 1, 0)

#N_sg_stima <- N_SUBGROUP (Deci_stima,         "tree_sg")

#discovery_stima <- DISCOVERY(N_sg_stima,  N_sg_true,   TF_D_stima)

#discovery_stima,




#pred_l2svm <- predict(l2svm)
#head30 <- as.data.frame( head(pred_l2svm$data, n=30) )
#tail30 <- as.data.frame( tail(pred_l2svm$data, n=30) )
#remove columns without constant value, if only 1 column keep its' column name
#head_con_var <- head30[,apply(head30, 2, var, na.rm=TRUE) == 0, drop=FALSE]
#tail_con_var <- tail30[,apply(tail30, 2, var, na.rm=TRUE) == 0, drop=FALSE]
#remove rownames (ID info)
#rownames(head_con_var)=c()
#rownames(tail_con_var)=c()
#when no HTE, InBoth="Treatment.effect", which is constant = -0.42
#InBoth = intersect(colnames(head_con_var), colnames(tail_con_var))
#when no HTE, HT_con_var= 30*2 observations of DeltaY (-0.42)
#HT_con_var =rbind(head_con_var[,InBoth], tail_con_var[,InBoth])
#when no HTE, INT_l2svm=NULL
#if ( is.null(colnames(HT_con_var))) {
#  INT_l2svm = "NA" } else {
#    INT_l2svm <- colnames(HT_con_var)
#  }





#mmatrix_rwl <- model.matrix(Y~(W + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)^4,Train_m)
#f <-reformulate( colnames(mmatrix_rwl)[-1] )

#f <- as.formula( ~(X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X1:X3) )


#f <- as.formula(paste("~", paste(names(Train_m)[!names(Train_m) %in% c("z.a_trueps", "W", "Y1", "Y", "fm")], collapse="+")))
#moPropen <- buildModelObj(model = ~ 1, solver.method = "glm", 
#                          solver.args = list('family'='binomial'),
#                          predict.method = 'predict.glm',
#                          predict.args = list(type = 'response'))
#moMain = buildModelObj(model = f, solver.method = "lm")
#txName = 'W'

#fitRWL <- rwl(moPropen = moPropen, moMain = moMain, 
#              data = Train_m %>% select(-c("z.a_trueps", "Y1", "fm")), 
#              reward = Train_m$Y,
#              txName = txName, regime = f, verbose = 0L, lambdas = 1)

##Available methods


#library(causalTree)
#Y_form <- as.formula(paste0("Y ~ ", paste0(colnames(X), collapse = " + ")))
#causaltree <- causalTree::causalTree(Y_form, data = Train_m, treatment = Train_m$W,
#                                     split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
#                                     xval = 5, cp = 0, 
#                                     minsize = 300, #requires at least minsize treated and minsize control observations in each leaf
#                                     propensity = 0.5)
#opcp <- causaltree$cptable[,1][which.min(causaltree$cptable[,4])]
#opfit <- prune(causaltree, opcp)
#rpart.plot(opfit)
#str(opfit)
#sp <- as.data.frame(opfit$splits) 
#splits <- sp %>% tibble::rownames_to_column() %>%
#  rename(n=count) %>%
#  mutate(var = ifelse( grepl("\\.", rowname)==T, gsub("\\..*","", rowname) , rowname ) )
#dplyr::left_join(opfit$frame, dplyr::select(splits,  c(var, n, ncat, improve, index, adj)), 
#                 by = c("var","n") )




# --------------------------
# STIMA
# --------------------------
#s_stima =  Sys.time()
#stima_fit<- stima::stima( Train_m, #be careful! 1st column is outcome, remaining as predictors!!!
#                          maxsplit=7, #each compute splitting node is time consuming, so only use 7
#                          first=3, #the column number in the data frame of the predictor that is used for the first split of the regression trunk. The default option automatically selects the predictor for the first split
#                          vfold=10)
#plot(stima_fit,digits=1)
#summary(stima_fit)
#prune the regression trunk
#stima_pr<- prune(stima_fit,data=Train_m)
#plot(stima_pr,digits=1)
#btree_stima <- stima_pr$trunk

#if( (is.null(btree_stima)==FALSE & nrow(btree_stima)==1) || is.null(btree_stima)==TRUE ) {
#  stima_subgroup="NA" 
#  stima_subgroup_key="NA" 
#  vote_stima_subgroup="NA" 
#} else {
#  stima_btree_df <- CONVERT_tree_stima2CF(btree_stima)
#  stima_subgroup     <- PRE_MAJORITY_SUBGROUP ( list(stima_btree_df) );  
#  stima_subgroup_key <- lapply(stima_subgroup, function(df) subset(df, select=c("subgroupID", "subgroup")))
#  vote_stima_subgroup<-MAJORITY_VOTE(list(stima_btree_df),  stima_subgroup_key,   stima_subgroup,   stima_subgroup, #similar as IT tree has not split frequency, this IT_subgroup is just to fool the Fx
#                                     tree_true_subgroup, tree_true_N1,   tree_true_N123)
#}
#e_stima =  Sys.time()
#time_stima = e_stima - s_stima