
#' DATA4grplasso
#' 
#' Function that prepare training or testing data for model selection based on D5, D4,D3,and D2 subgroup decisions 
#' @param V_D5_subgroup subgroup decision from depth 5 CF
#' @param V_D4_subgroup subgroup decision from depth 4 CF
#' @param V_D3_subgroup subgroup decision from depth 3 CF 
#' @param V_D2_subgroup subgroup decision from depth 2 CF 
#' @param dat the dataset 
#' 
#' @return the data and formula for group lasso
#' 
#' @export

DATA4grplasso <- function(V_D5_subgroup, V_D4_subgroup, V_D3_subgroup, V_D2_subgroup, dat){
 # Function that obtain subgroup ID 
 Dat_ID_SG_D5 <- GET_SUBGROUP_ID(V_D5_subgroup$majority, dat)
 Dat_ID_SG_D4 <- GET_SUBGROUP_ID(V_D4_subgroup$majority, dat)
 Dat_ID_SG_D3 <- GET_SUBGROUP_ID(V_D3_subgroup$majority, dat)
 Dat_ID_SG_D2 <- GET_SUBGROUP_ID(V_D2_subgroup$majority, dat)
 
   Dat_ID_SG <-Dat_ID_SG_D5 %>% 
      dplyr::rename(G5_define = Definition, 
                    G5        = SubgroupID)  %>%
      dplyr::mutate(G5_define = as.factor(G5_define),
                    G4        = Dat_ID_SG_D4$SubgroupID,
                    G4_define = as.factor(Dat_ID_SG_D4$Definition),
                    G3        = Dat_ID_SG_D3$SubgroupID,
                    G3_define = as.factor(Dat_ID_SG_D3$Definition),
                    G2        = Dat_ID_SG_D2$SubgroupID,
                    G2_define = as.factor(Dat_ID_SG_D2$Definition)) 
   
    formula_g2  <-  as.formula(  paste0("Y ~ W + G2 + W:G2 + ",                          paste0(colnames(X), collapse = " + ") )   )
    formula_g3  <-  as.formula(  paste0("Y ~ W + G3 + W:G3 + ",                          paste0(colnames(X), collapse = " + ") )   )
    formula_g4  <-  as.formula(  paste0("Y ~ W + G4 + W:G4 + ",                          paste0(colnames(X), collapse = " + ") )   )
    formula_g5  <-  as.formula(  paste0("Y ~ W + G5 + W:G5 + ",                          paste0(colnames(X), collapse = " + ") )   )

  return(list(Dat_ID_SG   = Dat_ID_SG, 
              formula_g2  = formula_g2,
              formula_g3  = formula_g3,
              formula_g4  = formula_g4,
              formula_g5  = formula_g5) )
}
