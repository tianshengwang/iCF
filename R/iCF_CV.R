

#' function to run Cross-Validation on whole dataset (except for external validation set) to get final subgroup decision G_iCF
#' @param dat dataset for overall population
#' @param K cross-validation fold 
#' @treeNo tree No
#' @iterationNo iteration No (if = 1 then oneCF)
#'  
#' @return final subgroup decision G_iCF
#' 
#' @export
#' 
#'

iCFCV <- function(dat, K, treeNo, iterationNo, min.split.var, split_val_round_posi, 
                  #tune_unit, 
                  P_threshold, variable_type, hdpct, HTE_P_cf.raw){
  
  s_iCF_CV=Sys.time()
  
  if ( round(HTE_P_cf.raw,1) <= P_threshold ){
    
    model.g2.act <- list()
    model.g3.act <- list()
    model.g4.act <- list()
    model.g5.act <- list()
    
    model.m.act  <- list()
    model.g2.tran <- list()
    model.g3.tran <- list()
    model.g4.tran <- list()
    model.g5.tran <- list()
    model.m.tran  <- list()
    
    cf_raw_key <- list()
    mse_per_fold.g2.act     <- list()
    mse_per_fold.g3.act     <- list()
    mse_per_fold.g4.act     <- list()
    mse_per_fold.g5.act     <- list()
    mse_per_fold.m.act      <- list()
    mse_per_fold.g2.tran    <- list()
    mse_per_fold.g3.tran    <- list()
    mse_per_fold.g4.tran    <- list()
    mse_per_fold.g5.tran    <- list()
    mse_per_fold.m.tran     <- list()
    
    SG_D <- list()
    vote_D5_tree.syn <- list()
    vote_D4_tree.syn <- list()
    vote_D3_tree.syn <- list()
    vote_D2_tree.syn <- list()
    
    vote_D5_subgroup.L <- list()
    vote_D4_subgroup.L <- list()
    vote_D3_subgroup.L <- list()
    vote_D2_subgroup.L <- list()
    
    stability_D5_T_r <- list()
    stability_D4_T_r <- list()
    stability_D3_T_r <- list()
    stability_D2_T_r <- list()
    
    Deci_Final_iCF.act.tr <- list()
    Deci_Final_iCF.tran.tr <- list()
    Deci_Final_iCF.act.te <- list()
    Deci_Final_iCF.tran.te <- list()
    
    test_data.tran <- list()

    Test_ID_SG_iCF <- list()
    
    selected_cf.idx.L <- list()
    accuracy_per_fold <- list()

    set.seed(1234)
    tt_indicies <- caret::createFolds(y=dat[,1], k= K)

    for(f in 1:length(tt_indicies)){
      Train_cf <- dat[-tt_indicies[[f]],]
      Test_cf  <- dat[tt_indicies[[f]],]
      ID_cf <-1:nrow(dat)
      Train_ID_cf <- cbind(Train_cf, as.vector(ID_cf[-tt_indicies[[f]]])) %>% dplyr::rename (ID=`as.vector(ID_cf[-tt_indicies[[f]]])`)
      Test_ID_cf  <- cbind(Test_cf,  as.vector(       tt_indicies[[f]]) ) %>% dplyr::rename (ID=`as.vector(tt_indicies[[f]])`)  
      #dplyr::all_equal(Train_ID_cf, Test_ID_cf)
      #============================== raw full CF for CV training data ==============================
      #for each training set, X, Y, W ... are regenerated
      cf_raw_key[[f]] <- CF_RAW_key(Train_cf, min.split.var, variable_type, hdpct)   

      selected_cf.idx.L[[f]]  <- cf_raw_key[[f]]$selected_cf.idx

      #==============================     iCF for CV training data   ==============================
      
      #check the top of this file: SUBGROUP_PIPLINE
      #to get subgroup decisions G_D, and models build from G_D to predict Y*
      SG_D[[f]]<-SUBGROUP_PIPELINE(X = cf_raw_key[[f]]$X,
                                  Y = cf_raw_key[[f]]$Y, 
                                  W = cf_raw_key[[f]]$W, 
                                  Y.hat =  cf_raw_key[[f]]$Y.hat, 
                                  W.hat =  cf_raw_key[[f]]$W.hat, 
                                  selected_cf.idx =selected_cf.idx.L[[f]],
                                  leafsize, 
                                  treeNo, 
                                  iterationNo, 
                                  Ntrain = nrow(Train_cf), 
                                  Train_ID = Train_ID_cf ,
                                  Test_ID =  Test_ID_cf,
                                  variable_type,
                                  HTE_P_cf.raw,
                                  P_threshold)
      
      vote_D5_tree.syn[[f]] <- SG_D[[f]]$vote_D5_tree.syn
      vote_D4_tree.syn[[f]] <- SG_D[[f]]$vote_D4_tree.syn
      vote_D3_tree.syn[[f]] <- SG_D[[f]]$vote_D3_tree.syn
      vote_D2_tree.syn[[f]] <- SG_D[[f]]$vote_D2_tree.syn
      
      vote_D5_subgroup.L[[f]] <- (SG_D[[f]]$vote_D5_subgroup)
      vote_D4_subgroup.L[[f]] <- (SG_D[[f]]$vote_D4_subgroup)
      vote_D3_subgroup.L[[f]] <- (SG_D[[f]]$vote_D3_subgroup)
      vote_D2_subgroup.L[[f]] <- (SG_D[[f]]$vote_D2_subgroup)
      
      stability_D5_T_r[[f]] <- SG_D[[f]]$stability_D5_T_r
      stability_D4_T_r[[f]] <- SG_D[[f]]$stability_D4_T_r
      stability_D3_T_r[[f]] <- SG_D[[f]]$stability_D3_T_r
      stability_D2_T_r[[f]] <- SG_D[[f]]$stability_D2_T_r
      
      Deci_Final_iCF.tran.tr[[f]] <- SG_D[[f]]$Deci_Final_iCF.tran.tr
      
      model.g2.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g2
      model.g3.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g3
      model.g4.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g4
      model.g5.tran[[f]] =  Deci_Final_iCF.tran.tr[[f]]$model.g5
      model.m.tran[[f]]  =  Deci_Final_iCF.tran.tr[[f]]$model.basic
      
      Test_ID_SG_iCF[[f]] <- SG_D[[f]]$Test_ID_SG_iCF
      
      #' Function that does prepare data that 1) ready to apply the formulas including interacting terms between treatment and factor subgroup deicision variable; 
      #' #2)transformed outcome
      #test_data.act[[f]] = SGMODEL_DATA(Test_ID_SG_iCF[[f]], "actual")$dat_ID_SG_df
      test_data.tran[[f]] =SGMODEL_DATA(Test_ID_SG_iCF[[f]], "transform")$dat_ID_SG_df
      
      
      #========================================
      Deci_Final_iCF.tran.te[[f]] <- SG_D[[f]]$Deci_Final_iCF.tran.te
      
      
 #     if (length(unique(Y)) >=8){
      
      #get the number of G5, G4, G3, G2 columns
        Ncol_g2345= ncol(test_data.tran[[f]] %>% dplyr::select((contains( c("G5", "G4", "G3", "G2") ) ) ) %>%
                                                dplyr::select(starts_with("G")) %>% #when hdiCF, some variable (code) names include "G5-2", e.g. dx3_outpt_G47
                                                dplyr::select(ends_with(c("5", "4", "3", "2") ))
                                               )
        
  #GET PREDICTED VALUES FROM DIFFERENT MODELS      
          
  #--------------------------------------------------------------------------- --------------------
  #using lm to predict transformed outcome (continuous):
          test_data.tran[[f]]$predict.g2.tran <- predict(model.g2.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g2345)])
          test_data.tran[[f]]$predict.g3.tran <- predict(model.g3.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g2345)])
          test_data.tran[[f]]$predict.g4.tran <- predict(model.g4.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g2345)])
          test_data.tran[[f]]$predict.g5.tran <- predict(model.g5.tran[[f]], newdata=test_data.tran[[f]][,1:(2 + ncol(X) + Ncol_g2345)])
          test_data.tran[[f]]$predict.m.tran  <- predict(model.m.tran[[f]],  newdata=test_data.tran[[f]][,1:(2 + ncol(X))])
        #calcualte performance across resamples

        #Y_star is renamed as Y:
        mse_per_fold.g2.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g2.tran, obs = test_data.tran[[f]]$Y)
        mse_per_fold.g3.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g3.tran, obs = test_data.tran[[f]]$Y)
        mse_per_fold.g4.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g4.tran, obs = test_data.tran[[f]]$Y)
        mse_per_fold.g5.tran[[f]] <- caret::postResample(pred = test_data.tran[[f]]$predict.g5.tran, obs = test_data.tran[[f]]$Y)
        mse_per_fold.m.tran[[f]]  <- caret::postResample(pred = test_data.tran[[f]]$predict.m.tran,  obs = test_data.tran[[f]]$Y)
 #--------------------------------------------------------------------------- --------------------  
        
    } #f loop over
    
    vote_D5_subgroup.L <<- vote_D5_subgroup.L
    vote_D4_subgroup.L <<- vote_D4_subgroup.L
    vote_D3_subgroup.L <<- vote_D3_subgroup.L
    vote_D2_subgroup.L <<- vote_D2_subgroup.L
    
    stability_D2_T_r <<- stability_D2_T_r
    stability_D3_T_r <<- stability_D3_T_r
    stability_D4_T_r <<- stability_D4_T_r
    stability_D5_T_r <<- stability_D5_T_r
    ##############################  
    # f iteration ends
    ##############################  
      # Bind together, add MSE

      #This function returns the mean MSE value across CV for main effect model, has nothing to do with HTE or SG decision!
      #mse_all_folds.m.act  <- CVBIAS_MAIN  (mse_per_fold.m.act)
      mse_all_folds.m.tran <- CVBIAS_MAIN  (mse_per_fold.m.tran)
    
    
    cv_sg_majority_D2 <- CV_SG_MAJORITY(vote_D2_subgroup.L, stability_D2_T_r, K)
    cv_sg_majority_D3 <- CV_SG_MAJORITY(vote_D3_subgroup.L, stability_D3_T_r, K)
    cv_sg_majority_D4 <- CV_SG_MAJORITY(vote_D4_subgroup.L, stability_D4_T_r, K)
    cv_sg_majority_D5 <- CV_SG_MAJORITY(vote_D5_subgroup.L, stability_D5_T_r, K)
    
    #the stability across CV for D2,D3,D4,and D5 SG
      stability_D2_CV <<- round( cv_sg_majority_D2$stability_CV, 3)
      stability_D3_CV <<- round( cv_sg_majority_D3$stability_CV, 3)
      stability_D4_CV <<- round( cv_sg_majority_D4$stability_CV, 3)
      stability_D5_CV <<- round( cv_sg_majority_D5$stability_CV, 3)
    
    #the stability of final SG in the CV (if plurity vote applied in CV due to different SG, then this is for the voted SG in CV )
    #e.g. if 3 of 5 CV D2 SG are the same, then the stability is the mean of the 3 stability scores of the 3 D2 SG
      stability_D2_vote <<- round(mean(unlist(stability_D2_T_r [c( cv_sg_majority_D2$selected_CV_idx) ] )), 3)
      stability_D3_vote <<- round(mean(unlist(stability_D3_T_r [c( cv_sg_majority_D3$selected_CV_idx) ] )), 3)
      stability_D4_vote <<- round(mean(unlist(stability_D4_T_r [c( cv_sg_majority_D4$selected_CV_idx) ] )), 3)
      stability_D5_vote <<- round(mean(unlist(stability_D5_T_r [c( cv_sg_majority_D5$selected_CV_idx) ] )), 3)
      
      stability_D2345_vote = cbind(stability_D2_vote, stability_D3_vote, stability_D4_vote, stability_D5_vote) %>% knitr::kable()
      
      stability_D2345_CV = cbind(stability_D2_CV, stability_D3_CV, stability_D4_CV, stability_D5_CV) %>% knitr::kable()
      
      stability_vote_across_CV = cbind(stability_D2_T_r, stability_D3_T_r, stability_D4_T_r, stability_D5_T_r) %>% as.data.frame() 
      #stability_vote_across_CV <- tibble::rowid_to_column(stability_vote_across_CV, "CVfold") 
      
      vote_D2_subgroup.L.majority <- cv_sg_majority_D2$majority_CV
      vote_D3_subgroup.L.majority <- cv_sg_majority_D3$majority_CV
      vote_D4_subgroup.L.majority <- cv_sg_majority_D4$majority_CV
      vote_D5_subgroup.L.majority <- cv_sg_majority_D5$majority_CV
    #This function returns the mean MSE value of those majority voted subgroup decisions across CV for depth=2,3,4,and 5

      mse_all_folds.g2.tran <- CVBIAS_D_MAJORITY  (vote_D2_subgroup.L, stability_D2_T_r, K, mse_per_fold.g2.tran)
      mse_all_folds.g3.tran <- CVBIAS_D_MAJORITY  (vote_D3_subgroup.L, stability_D3_T_r, K, mse_per_fold.g3.tran)
      mse_all_folds.g4.tran <- CVBIAS_D_MAJORITY  (vote_D4_subgroup.L, stability_D4_T_r, K, mse_per_fold.g4.tran)
      mse_all_folds.g5.tran <- CVBIAS_D_MAJORITY  (vote_D5_subgroup.L, stability_D5_T_r, K, mse_per_fold.g5.tran)
      
     #inverse CV stability score weight
      mse_all_folds.g2.tran_icvsw <- mse_all_folds.g2.tran/stability_D2_CV 
      mse_all_folds.g3.tran_icvsw <- mse_all_folds.g3.tran/stability_D3_CV
      mse_all_folds.g4.tran_icvsw <- mse_all_folds.g4.tran/stability_D4_CV
      mse_all_folds.g5.tran_icvsw <- mse_all_folds.g5.tran/stability_D5_CV
      


      ATE_all <- CATE_SG(dat, "NA")
      ATE_kable <- ATE_all %>% dplyr::mutate(SubgroupID="NA", Definition="Overall Population") #%>% knitr::kable()
      #--------------------------------------------------------------------------------------------------------------------------------------------------
      #Original, without any stability weight
      PICK_CV_ori <- PICK_CV(HTE_P_cf.raw, mse_all_folds.m.tran, 
                             mse_all_folds.g2.tran, mse_all_folds.g3.tran, mse_all_folds.g4.tran, mse_all_folds.g5.tran, 
                             vote_D2_subgroup.L, vote_D3_subgroup.L, vote_D4_subgroup.L, vote_D5_subgroup.L,
                             stability_D2_T_r, stability_D3_T_r,  stability_D4_T_r,  stability_D5_T_r,
                              P_threshold, K )
      CATE_ori <- CATE_SG(dat, PICK_CV_ori$selectedSG$majority)
      CATE_ori_kable <-  CATE_ori %>% mutate(Wt="NA") #%>% knitr::kable()
      #--------------------------------------------------------------------------------------------------------------------------------------------------
      #Inverse CV stability weight
      PICK_CV_icvsw <- PICK_CV(HTE_P_cf.raw, mse_all_folds.m.tran, 
                             mse_all_folds.g2.tran_icvsw, mse_all_folds.g3.tran_icvsw, mse_all_folds.g4.tran_icvsw, mse_all_folds.g5.tran_icvsw, 
                             vote_D2_subgroup.L, vote_D3_subgroup.L, vote_D4_subgroup.L, vote_D5_subgroup.L,
                             stability_D2_T_r, stability_D3_T_r,  stability_D4_T_r,  stability_D5_T_r,
                             P_threshold, K )
      
            if (identical(PICK_CV_icvsw$selectedSG,PICK_CV_ori$selectedSG)==T ) {
      CATE_icvsw <- CATE_ori 
      } else {
      CATE_icvsw <- CATE_SG(dat, PICK_CV_icvsw$selectedSG$majority)
      }
      CATE_icvsw_kable <-  CATE_icvsw %>% mutate(Wt="icvsw") #%>% knitr::kable()
 
    
  } else if (round(HTE_P_cf.raw,1) > P_threshold){
    
    time_iCFCV = "NA"
    
    ATE_all <- CATE_SG(dat, "NA")
    ATE_kable <- ATE_all %>% dplyr::mutate(SubgroupID="NA", Definition="Overall Population") %>% knitr::kable()
    
    vote_D2_subgroup.CVmajority = "NA"
    vote_D3_subgroup.CVmajority = "NA"
    vote_D4_subgroup.CVmajority = "NA"
    vote_D5_subgroup.CVmajority = "NA"
    
    stability_D2345_CV = "NA"
    
    stability_D2345_vote = "NA"
    
    stability_vote_across_CV = "NA"
    
    selectedSG_ori             = "NA"
    selectedSG_stability_ori   = "NA"
    CATE_t2_ori = CATE_ori_kable  = "NA"
    
    selectedSG_icvsw             = "NA"
    selectedSG_stability_icvsw    = "NA"
    CATE_t2_icvsw = CATE_icvsw_kable  = "NA"
    
    vote_D2_SG.L  = "NA"
    vote_D3_SG.L  = "NA"
    vote_D4_SG.L  = "NA"
    vote_D5_SG.L  = "NA"
    
  }
  
  e_iCF_CV=Sys.time()
  time_iCF_CV = difftime(e_iCF_CV, s_iCF_CV, units="secs")
  
  return(list(time_iCFCV = time_iCF_CV,
              
              ATE_kable = ATE_kable,
              
              vote_D2_subgroup.CVmajority = vote_D2_subgroup.L.majority,
              vote_D3_subgroup.CVmajority = vote_D3_subgroup.L.majority,
              vote_D4_subgroup.CVmajority = vote_D4_subgroup.L.majority,
              vote_D5_subgroup.CVmajority = vote_D5_subgroup.L.majority,
              
              stability_D2345_CV = stability_D2345_CV ,
              
              stability_D2345_vote = stability_D2345_vote ,
              
              stability_vote_across_CV = stability_vote_across_CV,
              
              selectedSG_ori             = PICK_CV_ori$selectedSG,
              selectedSG_stability_ori   = PICK_CV_ori$selectedSG_wt, 
              CATE_t2_ori = CATE_ori_kable ,
              
              selectedSG_icvsw             = PICK_CV_icvsw$selectedSG ,
              selectedSG_stability_icvsw   = PICK_CV_icvsw$selectedSG_wt, 
              CATE_t2_icvsw = CATE_icvsw_kable ,

              vote_D2_SG.L = vote_D2_subgroup.L %>% knitr::kable(),
              vote_D3_SG.L = vote_D3_subgroup.L %>% knitr::kable(),
              vote_D4_SG.L = vote_D4_subgroup.L %>% knitr::kable(),
              vote_D5_SG.L = vote_D5_subgroup.L %>% knitr::kable()
              

              )
         )
  
}


