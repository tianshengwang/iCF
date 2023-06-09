#-------------------------------
# Variable Importance
#-------------------------------
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


GG_VI <- function(varimp_cf, title){
  #VI_label <- colnames(X)
  varimp_cf_df <- as.data.frame( cbind( as.vector(varimp_cf), as.vector(colnames(X)) ) ) %>% 
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


#highlisght selected important variable in X-axis
#scale_x_discrete(labels=c("X1" =expression(bold(X1)), 
#"X3" =expression(bold(X3)),
#"X8"=expression(bold(X8)),
#"X2" =expression(bold(X2)), 
#"X4" =expression(bold(X4)),
#"X10"=expression(bold(X10)),
#parse=TRUE))

