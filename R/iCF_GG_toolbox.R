
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
  varimp_cf_df <- as.data.frame( cbind( as.vector(varimp_cf), as.vector( var_label) ) ) %>% 
    dplyr::rename(impValue=V1, vars_f=V2) %>%
    dplyr::mutate(impValue=as.character(impValue),
                  vars_f=as.character(vars_f)) %>%
    dplyr::mutate(impValue=as.numeric(impValue)) %>%
    dplyr::arrange(desc(impValue))
  
  VarImp <- varimp_cf_df %>%
    ggplot2::ggplot() +
    ggplot2::geom_bar(aes(y=reorder(vars_f, impValue),x=impValue, fill = impValue), stat = 'identity') + 
    ggplot2::scale_fill_viridis_c(name = "value") +
    ggplot2::theme_minimal() + 
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggplot2::labs(
      title = title,
      x = 'Importance Value',
      y= "Variables") +
    ggplot2::theme(plot.title    = element_text(color = "black", size=18, face="bold", #hjust = 0.5, 
                                       vjust=5, margin =ggplot2:: margin(0.5, 0, 0, 0, "cm")),
          plot.subtitle = element_text(color = "black", size=16, face="bold" ),
          axis.title.x  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
          axis.title.y  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
          plot.caption  = element_text(color = "black", size=12, face="plain"),
          axis.text.x   = element_text(color = "black", size=11, face="plain", angle = 0),
          axis.text.y   = element_text(color = "black", size=11, face="plain", angle = 0)) 
  
  
return(VarImp)
}


#' GG_PS
#' 
#' Function that shows Propensity Score distribution
#' @param dat dataset
#' @param PS predicted PS
#' @filename file name of the image 
#'  
#' @return the PS distribution
#' 
#' @export

GG_PS <- function(dat, PS, title, binwid, filename){
  dat_new = cbind(dat, PS)
  dat_new_W1 = subset(dat_new,W==1)
  dat_new_W0 = subset(dat_new,W==0)
png(filename=paste0(filename, ".png"))
PSplot <-ggplot2::ggplot(data=dat_new,
                      aes(x=PS, group=factor(W), fill=factor(W)))+
 # ggplot2:: geom_density(aes(y=..density..), alpha = 0.75, position = position_dodge(width=0.01)) +
  ggplot2::geom_histogram(aes(y=#..count.. 
                       ..density..
                     ),
                 alpha = 0.75,
                 binwidth=binwid,
                 position = position_dodge(width=0.01))+
  
  theme_classic()+
  xlab(paste(title, sep=": "))+ylab("Density")+
  labs(fill = "Treatment",
       title = paste0("Distribtution of ", title),
       subtitle = paste0("Number of variables: " , (ncol(dat)-2),
                         "; \n" ,
                         "Treatment=1: ",
                         "Min=" ,      round(summary( dat_new_W1$PS )[[1]] ,2 ),    
                         ";  Q1=",     round(summary( dat_new_W1$PS )[[2]] ,2 ),
                         ";  Median=", round(summary( dat_new_W1$PS )[[3]] ,2 ),
                         ";  Mean=",   round(summary( dat_new_W1$PS )[[4]] ,2 ),
                         ";  Q3=",     round(summary( dat_new_W1$PS )[[5]] ,2 ),
                         ";  Max=",    round(summary( dat_new_W1$PS )[[6]] ,2 ), 
                         "; \n" ,
                         "Treatment=0: ",
                         "Min=" ,      round(summary( dat_new_W0$PS )[[1]] ,2 ),    
                         ";  Q1=",     round(summary( dat_new_W0$PS )[[2]] ,2 ),
                         ";  Median=", round(summary( dat_new_W0$PS )[[3]] ,2 ),
                         ";  Mean=",   round(summary( dat_new_W0$PS )[[4]] ,2 ),
                         ";  Q3=",     round(summary( dat_new_W0$PS )[[5]] ,2 ),
                         ";  Max=",    round(summary( dat_new_W0$PS )[[6]] ,2 )
       ),
       caption = paste0("Sample size: ", length(PS) ) )

print(PSplot)
dev.off()

return(PSplot)
}


#' GG_IPTW
#' 
#' Function that shows IPTW distribution
#' @param dat dataset
#' @param IPTW IPTW
#'  
#' @return the PS distribution
#' 
#' @export

GG_IPTW <- function(IPTW){
IPTW_distri <- ggplot2::ggplot(data.frame(iptw_u= IPTW ), aes(iptw_u)) +   
               ggplot2::geom_histogram(aes(y=#..count.. 
                                             ..density..
                                           ),
                                           binwidth = 0.01,  fill='#62C6F2') +
               ggplot2::xlab(paste("Inverse probability weight (unstabilized)", sep=": "))+ylab("Density")+
               ggplot2::labs(title = "Distribution of Inverse Probability Weight",
                            subtitle = paste0("Min=" ,      round(summary( IPTW )[[1]] ,1 ),    
                                              ";  Q1=",     round(summary( IPTW )[[2]] ,1 ),
                                              ";  Median=", round(summary( IPTW )[[3]] ,1 ),
                                              ";  Mean=",   round(summary( IPTW )[[4]] ,1 ),
                                              ";  Q3=",     round(summary( IPTW )[[5]] ,1 ),
                                              ";  Max=",    round(summary( IPTW )[[6]] ,1 )),
                            
                            caption = paste0("Sample size: ", length(IPTW) ) )
                
return(IPTW_distri)
}




#' GG_Y_star
#' 
#' Function that shows transformed outcome distribution
#' @param dat dataset
#' @param IPTW IPTW
#'  
#' @return the PS distribution
#' 
#' @export

GG_Y_star <- function(dat, Y_star,No0 ){
  if (No0==F){
    No0flag=""
    dat_new=cbind(dat, Y_star)
  } else {
    No0flag="Non-zero "
    dat_new=cbind(subset(dat, Y==1), Y_star)
  }
  
  
  dat_new_W1 = subset(dat_new,W==1)
  dat_new_W0 = subset(dat_new,W==0)
  
  Y_star_distri <- ggplot2::ggplot( dat_new #data.frame(Transformed_O = Y_star ) 
                                   
                                   , aes(Y_star #Transformed_O
                                         ) ) +   
    ggplot2::geom_histogram(aes(y=#..count.. 
                                  ..density.., 
                                group=factor(W), fill=factor(W)
                                    ),
                            alpha = 0.7,
                            binwidth = 0.1,#,  fill='#62C6F2'
                            position = position_dodge(width=0.01)) +
    facet_grid(W ~ ., margins=F) +
    ggplot2::xlab(paste0(No0flag, "Transformed Outcome"))+ylab("Density")+
    ggplot2::labs(fill = "Treatment",
                  title = paste0("Distribution of ",No0flag, "Transformed Outcome"),
                  subtitle = paste0("Number of variables: " , (ncol(dat)-2),
                                    "; \n" ,
                                    "Treatment=1: ",
                                    "Min=" ,      round(summary(dat_new_W1$Y_star)[[1]] ,3 ),    
                                    ";  Q1=",     round(summary(dat_new_W1$Y_star)[[2]] ,3 ),
                                    ";  Median=", round(summary(dat_new_W1$Y_star)[[3]] ,3 ),
                                    ";  Mean=",   round(summary(dat_new_W1$Y_star)[[4]] ,3 ),
                                    ";  Q3=",     round(summary(dat_new_W1$Y_star)[[5]] ,3 ),
                                    ";  Max=",    round(summary(dat_new_W1$Y_star)[[6]] ,3 ),
                                    "; \n" ,
                                    "Treatment=0: ",
                                    "Min=" ,      round(summary(dat_new_W0$Y_star)[[1]] ,3 ),    
                                    ";  Q1=",     round(summary(dat_new_W0$Y_star)[[2]] ,3 ),
                                    ";  Median=", round(summary(dat_new_W0$Y_star)[[3]] ,3 ),
                                    ";  Mean=",   round(summary(dat_new_W0$Y_star)[[4]] ,3 ),
                                    ";  Q3=",     round(summary(dat_new_W0$Y_star)[[5]] ,3 ),
                                    ";  Max=",    round(summary(dat_new_W0$Y_star)[[6]] ,3 )) ,
                  caption =  paste0("Number of obserbations: ", length(Y_star))
                      )
  return(Y_star_distri)
}


#' ROC2comp
#' 
#' Function that compares 2 ROCs, not necessary for iCF/hdiCF
#' @param response1 true outcome 1
#' @param response2 true outcome 2
#' @param predict1 prediction 1
#' @param predict2 prediction 2  
#' @param label1 label 1
#' @param label2 lebel 2 
#' @param outcome outcome, e.g. Y or W
#' @return the image object
#' 
#' @export


ROC2comp <- function(response1, predict1, response2, predict2,label1, label2, outcome){
dev.off()

par(pty = "s") 
roc_Y.hat_all <-pROC::roc(response = response1, predictor = predict1, plot=T, 
                          legacy.axes=TRUE, percent=TRUE, 
                          xlab="False Positive Percentage", ylab="True Postive Percentage", 
                          col="#377eb8", lwd=4, print.auc=TRUE)

pROC::plot.roc(response2, predict2, percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", 
       legend=c(label1, label2), 
       col=c("#377eb8", "#4daf4a"), 
       lwd=4)
title(main =  paste0("ROC of predicted ", outcome) ,  line = 3, adj = 0)

ROC_compare <- recordPlot()
return( ROC_compare )

}