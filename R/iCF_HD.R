#' FIX_LOW_FREQ
#' data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==
#' Function that predict PS using LASSO algorithm, from the Personalized package,
#' need to modify as non HD variables may include more than 3 (age, sex, race)
#' @param dat matrix of features, X
#' @param Lcutoff cut off value for count as low frequnt (e.g. 10), lower than this value will lead to errors in cross validation 
#'  
#' @return a matrix of X removed low frequency (combined with the next lower level)
#' 
#' @export

FIX_LOW_FREQ <- function(dat, Lcutoff){
#-------------------------------------------------
# step 1. identify HD variables
#-------------------------------------------------
dat_HD <- dat %>% dplyr::select(starts_with(c( paste0("dx",dxgroup),"cpt5",paste0("atc", atcgroup) )))
  
##############################################################################
# step 2. reassign level value to the next level
##############################################################################
  
#-----------------------------------------------------------------------------
#  LEVEL 4 (X=3): assign those with extreme low freq value no the next level #
#-----------------------------------------------------------------------------
#kevel 4: 3 or frequent (> 75th percentile number of times).

X_3 <- lapply(dat_HD, function(x) (table(x) %>% as.data.frame())[4,2] )
min(unlist(X_3))
sort((unlist(X_3)))[1:25]
#level 4 with count < 10
HD_noLevel4 <- names((unlist(X_3))[unlist(X_3) < Lcutoff ])
#extract names of variables
HD_noLevel_4 <- HD_noLevel4[!is.na(HD_noLevel4)]
lapply(dat[HD_noLevel_4], function(x) table(x, useNA = "ifany") )


# for those with very low level 4 (value=3), when the HD variable=3, reassign to 2
for(var in HD_noLevel_4) {
  dat[[var]][dat[[var]] == 3] <- 2
}

lapply(dat[HD_noLevel_4], function(x) table(x, useNA = "ifany") )


#-----------------------------------------------------------------------------
#  LEVEL 3 (X=2): assign those with extreme low freq value no the next level #
#-----------------------------------------------------------------------------
#level 3: 2 sporadically (1 < appeared < 75th percentile number of times), 
X_2 <- lapply(dat_HD, function(x) (table(x) %>% as.data.frame())[3,2] )
min(unlist(X_2));
sort((unlist(X_2)))[1:15]
#level 3 with count < 10
HD_noLevel3 <- names((unlist(X_2))[unlist(X_2) < Lcutoff ])
#extract names of variables
HD_noLevel_3 <- HD_noLevel3[!is.na(HD_noLevel3)]
lapply(dat[HD_noLevel_3], function(x) table(x, useNA = "ifany") )


# for those missing level 3 (value=2), when the HD variable=3, reassign to 1
#after reassign values, level 4(X=3) may not exist, thus level 3 (X=2) may be the one with very low frequency
#those missing level 3 (X=2) will have only 3 levels: 0, 1, 3
for(var in HD_noLevel_3) {
  dat[[var]][dat[[var]] == 3] <- 1 
}

lapply(dat[HD_noLevel_3], function(x) table(x, useNA = "ifany") )
#------------------------------------------------------------------------------


#MUST use 2 separate steps, or cause errors (only one level remaining!!!)
#those after reassigning level 4 (X=3) to level 3 (X=2) that still very low may pop up!
for(var in HD_noLevel_3) {
  dat[[var]][dat[[var]] == 2] <- 1 
}

lapply(dat[HD_noLevel_3], function(x) table(x, useNA = "ifany") )


#-----------------------------------------------------------------------------
#  LEVEL 2: assign those with extreme low freq value no the next level #
#-----------------------------------------------------------------------------
#level 2: 1 occurred 1 time,  
X_1 <- lapply(dat_HD, function(x) (table(x) %>% as.data.frame())[2,2] )
min(unlist(X_1));
sort((unlist(X_1)))[1:15]
#level 2 with count < 10
HD_noLevel2 <- names((unlist(X_1))[unlist(X_1) < Lcutoff ])
#extract names of variables
HD_noLevel_2 <- HD_noLevel2[!is.na(HD_noLevel2)]
# for those with very low level 2 (value=1), when the HD variable=1, reassign to 0
lapply(dat[HD_noLevel_2], function(x) table(x, useNA = "ifany") )
#------------------------------------------------------------------------------

for(var in HD_noLevel_2) {
  dat[[var]][dat[[var]] == 1] <- 0
}

lapply(dat[HD_noLevel_2], function(x) table(x, useNA = "ifany") )

#-----------------------------------------------------------------------------
#  LEVEL 1: NO NEED TO WORRY #
#-----------------------------------------------------------------------------
X_0 <- lapply(dat_HD, function(x) (table(x) %>% as.data.frame())[1,2] )
min(unlist(X_0));
sort((unlist(X_0)))[1:15]


return(dat)
}

