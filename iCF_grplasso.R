
#' Function that prepare training or testing data for model selection based on D4,D3,and D2 subgroup decisions 
#' @param dat the dataset 
#' @param V_D4_subgroup subgroup decision from depth 4 CF
#' @param V_D3_subgroup subgroup decision from depth 3 CF 
#' @param V_D2_subgroup subgroup decision from depth 2 CF 
#' 
#' @return the data and formula for group lasso
#' 
#' @export


DATA4grplasso <- function(V_D4_subgroup, V_D3_subgroup, V_D2_subgroup, dat){
  #------------------------------------------
  #solution 1: LASSO for factor variables
  #------------------------------------------
 Dat_ID_SG_D4 <- GET_SUBGROUP_ID(V_D4_subgroup$majority, dat)
 Dat_ID_SG_D3 <- GET_SUBGROUP_ID(V_D3_subgroup$majority, dat)
 Dat_ID_SG_D2 <- GET_SUBGROUP_ID(V_D2_subgroup$majority, dat)
  #scenario 1: if D4_subgroup = D3_subgroup = D2_subgroup
  if( identical(V_D4_subgroup$majority, V_D3_subgroup$majority) &
      identical(V_D4_subgroup$majority, V_D2_subgroup$majority) ){
   Dat_ID_SG <-Dat_ID_SG_D2 %>% 
      dplyr::rename(G2_define = Definition, 
                    G2        = SubgroupID)  %>%
      dplyr::mutate(G2_define = as.factor(G2_define)) 
    
    formula_g2 <- formula_g3 <- formula_g4 <- formula_gl  <- formula_g23 <-  as.formula(  paste0("Y ~ W + G2 + W:G2 +", paste0(colnames(X), collapse = " + ") )   )
  } else if (
    #scenario 2: if D4_subgroup = D3_subgroup OR D4_subgroup = D2_subgroup
    identical(V_D4_subgroup$majority, V_D3_subgroup$majority) |  
    identical(V_D4_subgroup$majority, V_D2_subgroup$majority) ) {
   Dat_ID_SG <-Dat_ID_SG_D3 %>% 
      dplyr::rename(G3_define = Definition, 
                    G3        = SubgroupID)  %>%
      dplyr::mutate(G2        = Dat_ID_SG_D2$SubgroupID,
                    G2_define = as.factor(Dat_ID_SG_D2$Definition)) 
    formula_g23  <- formula_gl <-  as.formula(  paste0("Y ~ W + G2 + G3 +  W:G2 + W:G3 + ", paste0(colnames(X), collapse = " + ") )   )
                    formula_g2 <-  as.formula(  paste0("Y ~ W + G2 + W:G2 + ",              paste0(colnames(X), collapse = " + ") )   )
                    formula_g3 <-  as.formula(  paste0("Y ~ W + G3 + W:G3 + ",              paste0(colnames(X), collapse = " + ") )   )
                    formula_g4 <-  as.formula(  paste0("Y ~ W + G4 + W:G4 + ",              paste0(colnames(X), collapse = " + ") )   )
    
    #scenario 3: if D3_subgroup = D2_subgroup
  } else if ( identical(V_D3_subgroup$majority, V_D2_subgroup$majority) ) {
   Dat_ID_SG <-Dat_ID_SG_D4 %>% 
      dplyr::rename(G4_define = Definition, 
                    G4        = SubgroupID)  %>%
      dplyr::mutate(G2        = Dat_ID_SG_D2$SubgroupID,
                    G2_define = as.factor(Dat_ID_SG_D2$Definition)) 
    formula_gl  <-                                as.formula(  paste0("Y ~ W + G2 + G4 + W:G2 + W:G4 + ", paste0(colnames(X), collapse = " + ") )   )
    formula_g23 <- formula_g3  <- formula_g2  <-  as.formula(  paste0("Y ~ W + G2 + W:G2 + ",             paste0(colnames(X), collapse = " + ") )   )
    formula_g4  <-                                as.formula(  paste0("Y ~ W + G4 + W:G4 + ",             paste0(colnames(X), collapse = " + ") )   )
  } else {
    #scenario 4: none are the same
   Dat_ID_SG <-Dat_ID_SG_D4 %>% 
      dplyr::rename(G4_define = Definition, 
                    G4        = SubgroupID)  %>%
      dplyr::mutate(G4_define = as.factor(G4_define),
                    G3        = Dat_ID_SG_D3$SubgroupID,
                    G3_define = as.factor(Dat_ID_SG_D3$Definition),
                    G2        = Dat_ID_SG_D2$SubgroupID,
                    G2_define = as.factor(Dat_ID_SG_D2$Definition)) 
    formula_gl  <-  as.formula(  paste0("Y ~ W + G2 + G3 + G4 +  W:G2 + W:G3 + W:G4 + ", paste0(colnames(X), collapse = " + ") )   )
    formula_g23 <-  as.formula(  paste0("Y ~ W + G2 + G3 + W:G2 + W:G3  + ",             paste0(colnames(X), collapse = " + ") )   )
    formula_g2  <-  as.formula(  paste0("Y ~ W + G2 + W:G2 + ",                          paste0(colnames(X), collapse = " + ") )   )
    formula_g3  <-  as.formula(  paste0("Y ~ W + G3 + W:G3 + ",                          paste0(colnames(X), collapse = " + ") )   )
    formula_g4  <-  as.formula(  paste0("Y ~ W + G4 + W:G4 + ",                          paste0(colnames(X), collapse = " + ") )   )
  }
  
  return(list(Dat_ID_SG   = Dat_ID_SG, 
              formula_gl  = formula_gl,
              formula_g2  = formula_g2,
              formula_g3  = formula_g3,
              formula_g4  = formula_g4,
              formula_g23 = formula_g23) )
}




  
  
  
  
#' Function that run cross-validation for group lasso (grplasso) using a formula (formula_grplasso) designed to CF, 
#' BE CAREFUL!!! Need to revise for classification.
#' @param dat the dataset 
#' @param lambda_grplasso a grid of penalty parameter from 0 to lambda_max
#' 
#' @return the list of cross-validation results and the lambda leading to minimor. 
#' 
#' @export

#dat=Train_ID_SG_df 
#lambda_grplasso=lambda_grplasso_50
#formula_grplasso = formula_gl
#contrast = contr


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



#else if ( VS_method=="grplasso_g" ){
  #------------------------------------------
  #solution 2: GROUP LASSO for factor variables
  #------------------------------------------
# Dat_ID_SG_D4 <- GET_SUBGROUP_ID(vote_D4_subgroup$majority, Train_ID)
# Dat_ID_SG_D3 <- GET_SUBGROUP_ID(vote_D3_subgroup$majority, Train_ID)
# Dat_ID_SG_D2 <- GET_SUBGROUP_ID(vote_D2_subgroup$majority, Train_ID)
  
# Dat_ID_SG_D4 <- MAKE_GROUP_COLUMN(Train_ID_SG_D4, "D4",  "group", "g")
# Dat_ID_SG_D3 <- MAKE_GROUP_COLUMN(Train_ID_SG_D3, "D3",  "group", "g")
# Dat_ID_SG_D2 <- MAKE_GROUP_COLUMN(Train_ID_SG_D2, "D2",  "group", "g")
  #scenario 1: if D4_subgroup = D3_subgroup = D2_subgroup
#  if( identical(vote_D4_subgroup$majority, vote_D3_subgroup$majority) & identical(vote_D4_subgroup$majority, vote_D2_subgroup$majority) ){
#   Dat_ID_SG <-Dat_ID_SG_D2 %>% 
#     dplyr::rename(G_D2_define = Definition_D2, 
#                    G_D2 = SubgroupID_D2)  %>%
#      dplyr::mutate(G_D2_define = as.factor(G_D2_define))
#  } else if (
    #scenario 2: if D4_subgroup = D3_subgroup OR D4_subgroup = D2_subgroup
#    identical(vote_D4_subgroup$majority, vote_D3_subgroup$majority) |  identical(vote_D4_subgroup$majority, vote_D2_subgroup$majority) ) {
#   Dat_ID_SG <-Dat_ID_SG_D3 %>%  dplyr::left_join(Train_ID_SG_D2 %>% dplyr::select(-c("Y", "W", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")), by = "ID"  ) %>%
#      dplyr::rename(G_D3_define = Definition_D3, 
#                    G_D3 = SubgroupID_D3)  %>%
#      dplyr::mutate(G_D2 =          Dat_ID_SG_D2$SubgroupID_D2,
#                    G_D2_define = as.factor(Train_ID_SG_D2$Definition_D2)) 
    #scenario 3: if D3_subgroup = D2_subgroup
#  } else if ( identical(vote_D3_subgroup$majority, vote_D2_subgroup$majority) ) {
#   Dat_ID_SG <-Dat_ID_SG_D4 %>%  dplyr::left_join(Train_ID_SG_D2 %>% dplyr::select(-c("Y", "W", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")), by = "ID"  ) %>%
#      
#      dplyr::rename(G_D4_define = Definition_D4, 
#                    G_D4 = SubgroupID_D4)  %>%
#      dplyr::mutate(G_D2 =          Dat_ID_SG_D2$SubgroupID_D2,
#                    G_D2_define = as.factor(Train_ID_SG_D2$Definition_D2))
#  } else {
    #scenario 4: none are the same
#   Dat_ID_SG <-Dat_ID_SG_D4 %>%  dplyr::left_join(Train_ID_SG_D3 %>% dplyr::select(-c("Y", "W", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")), by = "ID"  ) %>%
#      dplyr::left_join(Train_ID_SG_D2 %>% dplyr::select(-c("Y", "W", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")), by = "ID"  ) %>%
#      dplyr::rename(G_D4_define = Definition_D4, 
#                    G_D4 = SubgroupID_D4)  %>%
#      dplyr::mutate(G_D3 =          Dat_ID_SG_D3$SubgroupID_D3,
#                    G_D3_define = as.factor(Train_ID_SG_D3$Definition_D3),
#                    G_D2 =          Dat_ID_SG_D2$SubgroupID_D2,
#                    G_D2_define = as.factor(Train_ID_SG_D2$Definition_D2))  
#  }  
#}