
#' GG_DepthDistribution
#' 
#' Function that merge distribution of depth for proposed denominator for minimumum leaf size for D2, D3, D4, and D5 
#' @param D2_distribtion distribution plot for D2
#' @param D3_distribtion distribution plot for D3
#' @param D4_distribtion distribution plot for D4
#' @param D5_distribtion distribution plot for D5
#'  
#' @return the prediced propensity score
#' 
#' @export


GG_DepthDistribution <- function(D2_distribtion,D3_distribtion,D4_distribtion, D5_distribtion, Name){
  tiff(paste(Name, ncol(X), ".tiff", sep="_"), units="in",width=15, height=25, res=150)
  g_depth <- cowplot::plot_grid( D2_distribtion,D3_distribtion,D4_distribtion, D5_distribtion, 
                                 ncol  = 2, nrow=2,
                                 labels = "AUTO", label_size = 20
  )
  dev.off()
  return(g_depth)
}


#-------------------------------
# Variable Importance
#-------------------------------
#' PlotVI
#' 
#' Function that make plot for variable imporce, add label and title
#' @param varimp_CF variable
#' @param title title text
#'  
#' @return the variable importance plot
#' 
#' @export

PlotVI <- function(varimp_CF,title){
 # VI_label <- colnames(X)
  plot(varimp_CF, main=paste0("Variable Importance for ", title),
       xlab="Index", ylab="Variable importance value",
  )
  text(1:ncol(X), 
       varimp_CF,
       labels=#colnames(X)
         colnames(X)
  )
}

#' GG_VI
#' 
#' Function that make bar plot for variable imporce, add label and title
#' @param varimp_CF variable
#' @param title title text
#' @param var_label title text
#'  
#' @return the variable importance plot
#' 
#' @export


GG_VI <- function(varimp_cf, title, var_label){
  #VI_label <- colnames(X)
  varimp_cf_df <- as.data.frame( cbind( as.vector(varimp_cf), as.vector( var_label
                                                                       # colnames(X)
                                                                         ) ) ) %>% 
    dplyr::rename(impValue=V1, vars_f=V2) %>%
    dplyr::mutate(impValue=as.character(impValue),
                  vars_f=as.character(vars_f)) %>%
    dplyr::mutate(impValue=as.numeric(impValue)) %>%
    dplyr::arrange(desc(impValue))
  
  VarImp <- varimp_cf_df %>%
    ggplot() +
    geom_bar(aes(y=reorder(vars_f, impValue),x=impValue, fill = impValue), stat = 'identity') + 
    scale_fill_viridis_c(name = "value") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(
      title = title,
      x = 'Importance Value',
      y= "Variables") +
    theme(plot.title    = element_text(color = "black", size=18, face="bold", #hjust = 0.5, 
                                       vjust=5, margin =ggplot2:: margin(0.5, 0, 0, 0, "cm")),
          plot.subtitle = element_text(color = "black", size=16, face="bold" ),
          axis.title.x  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
          axis.title.y  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
          plot.caption  = element_text(color = "black", size=12, face="plain"),
          axis.text.x   = element_text(color = "black", size=11, face="plain", angle = 0),
          axis.text.y   = element_text(color = "black", size=11, face="plain", angle = 0)) 
  
  
return(VarImp)
}



