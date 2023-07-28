rm(list=ls())

devtools::install_github("tianshengwang/iCF")

library(iCF)
#CF_RAW_key
#source("/local/iCF/best_tree_MSegar.R")
#PlotVI
#GG_VI
#MinLeafSizeTune
#iCF
#CF
#oneCF
#GET_TREE_DF
#GET_TREE_L 

source("/local/iCF/GG_toolbox.R")
source("/local/iCF/sim_GenSimData.R")
source("/local/iCF/best_tree_MSegar.R")
source("/local/iCF/iCF_TREE_build.R")
source("/local/iCF/iCF_PARENT_node.R")
source("/local/iCF/iCF_PRE_majority.R")
source("/local/iCF/sim_Truth_tree.R")
source("/local/iCF/iCF_MAJORITY_VOTE.R")
source("/local/iCF/iCF_SUBGROUP_DECISION.R")
source("/local/iCF/iCF_SG_PIPELINE.R")
source("/local/iCF/iCF_CV.R")
source("/local/iCF/iCF_SUBGROUP_ANALYSIS.R")


nstudy = 10000
b0 <- 0
b1 <- 0.8
b2 <- -0.25
b3 <- 0.6
b4 <- -0.4
b5 <- -0.8
b6 <- -0.5
b7 <- 0.7
a0 <- -3.85
a1 <- 0.3
a2 <- -0.36
a3 <- -0.73
a4 <- -0.2
a5 <- 0.71
a6 <- -0.19
a7 <- 0.26
g1 <- -0.4 

F.sample.cor <- function(X, rho) {
  Y <- (rho * (X - mean(X)))/sqrt(var(X)) + sqrt(1 - rho^2) * rnorm(length(X))
  return(Y)
}

X1 <- rnorm(nstudy, mean=0, sd=1)
X2 <- rnorm(nstudy, mean=0, sd=1)
X3 <- rnorm(nstudy, mean=0, sd=1)
X4 <- rnorm(nstudy, mean=0, sd=1)
X5 <- F.sample.cor(X1, 0.2)
X6 <- F.sample.cor(X2, 0.9)
X7 <- rnorm(nstudy, mean=0, sd=1)
X8 <- F.sample.cor(X3, 0.2)
X9 <- F.sample.cor(X4, 0.9)
X10 <- rnorm(nstudy, mean=0, sd=1)

X1 <- ifelse(X1 > mean(X1), 1, 0)
X3 <- ifelse(X3 > mean(X3), 1, 0)
X5 <- ifelse(X5 > mean(X5), 1, 0)
X6 <- ifelse(X6 > mean(X6), 1, 0)
X8 <- ifelse(X8 > mean(X8), 1, 0)
X9 <- ifelse(X9 > mean(X9), 1, 0)

PS <- (1 + exp( -(b0 + b1*X1 + b2*X2 + b3*X3 + b4*X4 + b5*X5 + b6*X6 + b7*X7) ))^-1 #true propensity score

W <- rbinom(nstudy,1,PS) #treatment assignment

#Y = a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*W + rnorm(nstudy,0,1) + 0.4*W*X3 
Y = a0 + a1*X1 + a2*X2 + a3*X3 + a4*X4 +a5*X8 + a6*X9 + a7*X10 + g1*W + rnorm(nstudy,0,1) + 0.4*W*X3 + 0.3*W*X1 + 0.4*W*X1*X3 + 0.2*X1*X3 
dat <<- as.data.frame(cbind(W, Y, X1, X2, X3 ,X4, X5, X6, X7, X8, X9, X10))



vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W" ))  )
intTRUE <- "Unknown"
X <- dat[,vars_forest]
Y <- as.vector( as.numeric( dat[,"Y"] ) )
W <- as.vector( as.numeric( dat[,"W"] ) )

cf_raw_key.tr <- CF_RAW_key(dat, 1, "non-hd", hdpct=0.90) 
Y.hat  <<- cf_raw_key.tr$Y.hat                 
W.hat  <<- cf_raw_key.tr$W.hat                 
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw    
varimp_cf  <- cf_raw_key.tr$varimp_cf          
length(W.hat); 
selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx 

GG_VI(varimp_cf, "Variable importance" )






split_val_round_posi=0
#Get truth
truth.list <<- TRUTH("Unknown")
##define variables categorical variables more than 2 levels
vars_catover2 <<- NA



D2_MLS=MinLeafSizeTune(dat=dat,denominator=25, treeNo = 1000, iterationNo=100, split_val_round_posi=0, depth= "D2")
D2_MLS$depth_mean
D2_MLS$depth_gg

D3_MLS=MinLeafSizeTune(dat=dat,denominator=45, treeNo = 1000, iterationNo=100, split_val_round_posi=0,  depth="D3")
D3_MLS$depth_mean
D3_MLS$depth_gg

D4_MLS=MinLeafSizeTune(dat=dat, denominator=65, treeNo = 1000, iterationNo=100, split_val_round_posi=0,  depth="D4")
D4_MLS$depth_mean
D4_MLS$depth_gg

D5_MLS=MinLeafSizeTune(dat=dat,denominator=85, treeNo = 1000, iterationNo=100, split_val_round_posi=0, depth="D5")
D5_MLS$depth_mean
D5_MLS$depth_gg



leafsize <<- list(D5=85, D4=65, D3=45, D2=25)

iCFCV_B1000_i200_sim <- iCFCV(dat=dat,
                              K=5,
                              treeNo=100, 
                              iterationNo=10,
                              min.split.var=4, 
                              split_val_round_posi=0, 
                              P_threshold=0.1, 
                              variable_type = "non-HD",
                              hdpct= 0.95,
                              HTE_P_cf.raw = HTE_P_cf.raw) 

iCFCV_B1000_i200_sim






setwd("/local/projects/medicare/DPP4i_HTE/data/R")
load("HHF_SGLTvGLP_iCFCV.RData")


dat00 <-  hfp_2yr_all_sgltvglp %>% 
          dplyr::filter(FillDate2 <= 21184 -365*2 &  #31DEC2017
                        IndexDate >= 19449 #10359# & 1/2/2013    
                        ) %>%
          dplyr::mutate(Y = hfp_2yr_2yr,
                        age = as.numeric( cut(age, c(65,70,75,80,85,Inf) ,
                                              labels=c("65<age<=70 ","70<age<=75","75<age<=80","80<age<=85", "age>85")
                        )
                        )
          ) %>% #defined in the beginning of analysis_MACE SAS macro
  dplyr::rename(W=SGLT) %>%
  dplyr::select(BENE_ID, IndexDate, FillDate2, Y, W,
                age, race2, 
                sex, bl_DIABRETINOPATHY,  bl_NEPHROPATHY,  bl_NEUROPATHY,   bl_DIABCOMPLICATIONS, bl_nodrugs,  bl_HYPERGLYCEMIA, 
                baselinecvd, baselinechf,
                bl_HYPOGLYCEMIA, bl_FOOTULCER, bl_angina2, bl_mi2, bl_INTERVENTION, bl_atherosclerosis2, bl_ischemichtdz2, bl_cerebrovasculardz2, 
                bl_cardiomyopathy2, bl_chf2,  bl_pvdnew, bl_AF, bl_ARRHYTHMIA, bl_CARDIACARREST, bl_DEFIBRILLATOR,
                bl_ANEMIA, bl_ALCOHOL,bl_ASTHMA, bl_BRAININJ, bl_CANCER, bl_CHRLUNG,bl_CKD, bl_COAGULOPATHY, bl_CONNECTIVE, bl_DEMENTIA, bl_DEFANEMIA, 
                bl_DEPRESSION, bl_DIFFWALK, bl_DYSLIPIDEMIA, bl_hypothyroidism, bl_EDEMA, bl_ELECTROLYTES, bl_HIV, bl_HEMATOLOGICAL, 
                bl_HYPERTENSION, bl_HYPOTENSION, bl_IMMUNE, bl_METABOLIC, bl_METASTATICCA, bl_MILDLIVER, bl_MODLIVER, bl_NUTRITIONAL, 
                bl_NERVOUS, bl_PARAPLEGIA, bl_PARKINSON, bl_PNEUMONIA,bl_PSYCHOSIS, bl_PULCIRC, bl_REHABILITATION, bl_RENAL, bl_RHEUMATIC, bl_smokingfinal, 
                bl_THROMBOEMBOLISM, bl_VULVULAR, bl_WTLOSS, bl_AMBULANCE, bl_HOSPITALBED, bl_OXYGEN, bl_WHEELCHAIR,  
                bl_METFORMIN, bl_SAINSULIN, bl_LAI, bl_TZD, bl_MEGLITINIDE, bl_SULF, bl_DPP, 
                bl_IMMUNOSUPPRESSIVE, bl_STEROIDS, bl_ACEI, bl_ARB, bl_CCB, bl_BB, bl_NSAIDS, bl_ASPIRIN, 
                bl_ORALCONTRAC, bl_ESTROGEN, bl_LOOP, bl_OTHERDIURETICS, bl_STATIN,
                bl_HBA1C,  bl_FLUSHOT , bl_LIPIDTEST,  bl_HOSP, bl_HOSPDAYS,bl_ED,bl_ERDM,bl_outpt, bl_outptdm
              ) %>% 
              as.data.frame.matrix() %>%
              mutate(sex=as.numeric(sex)) #%>%

vars_catover2 <- c( "race2",  "age", "bl_HOSP", "bl_HOSPDAYS", "bl_ED",   "bl_ERDM",    "bl_outpt",    "bl_outptdm")

truth.list <<- TRUTH("Unknown")

dat <<- dat00%>%  dplyr::select(-c("BENE_ID", "IndexDate", "FillDate2"))

Ntrain <- nrow(dat)
nrow(dat)
ID <-1:nrow(dat)


vars_forest = colnames( dat %>% dplyr::select(-c("Y", "W"))  ) 

X <<- dat[,vars_forest]
Y <<- as.vector( as.numeric( dat[,"Y"] ) )
W <<- as.vector( as.numeric( dat[,"W"] ) )


cf_raw_key.tr <- CF_RAW_key(dat, min.split.var=4,  variable_type="non-hd", hdpct=0.95)    

Y.hat  <<- cf_raw_key.tr$Y.hat
W.hat  <<- cf_raw_key.tr$W.hat
HTE_P_cf.raw <<- cf_raw_key.tr$HTE_P_cf.raw
varimp_cf  <- cf_raw_key.tr$varimp_cf
selected_cf.idx <<- cf_raw_key.tr$selected_cf.idx
time_rawCF <- cf_raw_key.tr$time_rawCF

X[,selected_cf.idx]


VI_label <- c( "age"    ,               "race"    ,             "sex"      ,            "diabetes retinopathy",       "diabetes nephropathy",
               "diabetes neuropathy",      "Diabetes circulatory complications", 
                "N of GLDs", "N of hyperglycemia diagnoses",
                "CVD" ,          "CHF",           "Hypoglycemia" ,      "Foot ulcers" ,         "Angina" ,          
                "MI",                "Cardiac revascularization or bypass",       "Atherosclerosis",   "Ischemic heart diseases",      "Cerebrovascular diseases", "Cardiomyopathy",   
                "Congestive heart failure",               "PVD",             "AF" ,                "Arrhythmia",         "Cardiac arrest",      "Defibrillator",     
                "Anemia",             "Alcohol disorders",            "Asthma",             "Brain injury",           "Cancer",             "Chronic lung disorders",           
                "CKD",                "Coagulopathy",       "Connective tissue disorders",         "Dementia",           "Deficiency anemia",          "Depression",        
                "Difficulty walking",           "Dyslipidemia",       "Hypothyroidism",     "Edema",              "Electrolytes disorders",       "HIV",               
                "Hematological disorders",      "Hypertension",       "Hypotension",        "Immune disorders",             "Metabolic disorders",          "Metastatic cancers",      
                "Mild liver disorders",          "Moderate liver disorders",           "Nutritional disorders",        "Nervous system disorders",            "Paraplegia",         "Parkinsonism",         
                "Pneumonia",          "Psychosis",          "Pulmonary circulation disorders",            "Rehabilitation",     "Renal disorders",              "Rheumatic disorders",         
                "Smoking and smoking cessation",       "Thromboembolism",    "Valvular disorders",           "Weight loss",             "Ambulance",          "Hospital beds",       
                "Home oxygen",             "Wheelchairs",         "Metformin",          "Short-acting insulin",          "LAI",                "TZD",               
                "Meglitinide",        "Sulfonylurea",               "DPP4i",                "Immunosuppressive drugs",  "Steroids",           "ACEI",              
                "ARB",                "CCB",                "BB",                 "NSAIDS",             "Aspirin",            "Oral contraceptives",       
                "Estrogen",           "LOOP",               "Other diuretics",     "Statin" ,
                "N of HbA1C tests", "N of flu shots", "N of lipid tests", "N of hospital admissions", "Days of hospitalization", 
                "N of emergency room visits","N of emergency room visits due to DM", 
                "N of outpatient visits", "Outpatient visits due to DM" 
              )


GG_VI(varimp_cf, 'Variable Importance for SGLT2i vs GLP1RA cohort for HFF', VI_label )

GG_VI(varimp_cf, 'Variable Importance for SGLT2i vs GLP1RA cohort for HFF', colnames(X) )



D2_MLS=MinLeafSizeTune(dat=dat, denominator=25, treeNo = 1000, iterationNo=100, split_val_round_posi=0, depth= "D2", "#62C6F2")
D2_MLS$depth_mean
D2_MLS$depth_gg

D3_MLS=MinLeafSizeTune(dat=dat, denominator=35, treeNo = 1000, iterationNo=100, split_val_round_posi=0,  depth="D3", "#62C6F2")
D3_MLS$depth_mean
D3_MLS$depth_gg

D4_MLS=MinLeafSizeTune(dat=dat, denominator=50, treeNo = 1000, iterationNo=100, split_val_round_posi=0,  depth="D4", "#62C6F2")
D4_MLS$depth_mean
D4_MLS$depth_gg

D5_MLS=MinLeafSizeTune(dat=dat, denominator=70, treeNo = 1000, iterationNo=100, split_val_round_posi=0, depth="D5","#62C6F2")
D5_MLS$depth_mean
D5_MLS$depth_gg


leafsize <<- list(D5=175, D4=125, D3=100, D2=60)
iCFCV_SG <- iCFCV(K=5,
                  treeNo=100, 
                  iterationNo=50,
                  min.split.var=4, 
                  P_threshold=0.1, 
                  variable_type = "non-hd",
                  hdpct=0.95, 
                  HTE_P_cf.raw = HTE_P_cf.raw
)

