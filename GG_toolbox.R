library(gridExtra)
library(viridis)
library(ggpubr)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggridges)
theme_set(theme_ridges())
library(stringr)
library(data.table)
library(purrr)
#------------------------------------


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

REMOVE_LEGEND <- function(input_gg) {
output_gg    <- 	input_gg     + theme(legend.position="none")
return(output_gg)
}


GG_MERGE_2h <- function(fig1,fig2,Name){
  tiff(paste(Name, treeNo, iterationNo, Ntrain, intTRUE, ".tiff", sep="_"), units="in",width=15, height=10, res=150)
  g4 <- cowplot::plot_grid( fig1, fig2,
                            ncol  = 2, labels = "AUTO", label_size = 20
  )
  dev.off()
  return(g4)
}


#---------------------------------------------------------------
#  SPLIT FREQUENCY HEAT MAP
#---------------------------------------------------------------



GG_SPLITFREQ_HEAT <- function(iCF_D_SF_df, tree_depth, type_i, iflegend, truth_description, vars, iterationNo){
  
  iCF_D_SF_df_NOK <- lapply(iCF_D_SF_df, function(df) df[ , !(names(df) %in% c("k"))])
  #distinguish between iteration=1, ALL, and majority iterations
  Stability <-  format(round(length(iCF_D_SF_df)/iterationNo, digits=2), nsmall = 2)    
#-----------------------------------------
# to distinguish raw, iCF, and voted iCF
#-----------------------------------------
    if (type_i=="i1"){
    #only 1 dataframe, no iteration!!!
    D = 	iCF_D_SF_df_NOK
    iteration_char = ""
    CF = "CF"
#    subTITLE = paste0(#"No. of forest: ", length(iCF_D_SF_df) 
#                      "Split Frequency"
#                      ) 
  } else if(type_i=="all") {
    #get mean of MATRIX, c(2,3) means transforms the data to a 3-dimensional array and then takes a column mean out of it
    D = plyr::aaply(plyr::laply(iCF_D_SF_df_NOK, as.matrix), c(2, 3), mean)	
    iteration_char = ""
    CF = "iCF"
#    subTITLE = paste0("No. of forests: ", length(iCF_D_SF_df) ) 
  } else if(type_i=="majority" & length(iCF_D_SF_df) == 1){
    D = 	iCF_D_SF_df_NOK
    iteration_char = "from voted best tree"
    CF= "iCF"
#    subTITLE = paste0("No. of forest: ", length(iCF_D_SF_df), "; stability score = ",  Stability)  
  } else {
    D = plyr::aaply(plyr::laply(iCF_D_SF_df_NOK, as.matrix), c(2, 3), mean)	
    iteration_char = "from voted best tree"
    CF="iCF"
#    subTITLE = paste0("No. of forests: ", length(iCF_D_SF_df), "; stability score = ",  Stability )  
  }

  if(tree_depth=="D2"){
    breaks_para <- c(-2,-1)
    title_para = ""
    subtitle_para = paste(tree_depth, CF, iteration_char)
    size_para=5
    label_para= c("2","1")
    
  } else if (tree_depth=="D3") {
    breaks_para <- c(-3,-2,-1)
    title_para = ""
    subtitle_para = paste(tree_depth, CF, iteration_char)
    size_para=5
    label_para= c("3","2","1")
    
    
  } else if (tree_depth=="D4") {
    breaks_para <- c(-4,-3, -2, -1)
    title_para = ""
    subtitle_para = paste(tree_depth, CF, iteration_char)
    size_para=5
    label_para= c("4", "3","2","1")
    
  } else {
    breaks_para <- c(-5, -4,-3, -2, -1)
    title_para = paste(truth_description)
    subtitle_para = "raw CF"
    size_para=4
    label_para= c("5", "4", "3","2","1")
    
  }
  
  if(iflegend=="nolegend"){
    legend_para = "none"} else {
    legend_para = "right"}
 

  d <- data.frame(D)
  p <- length(vars)
  max.depth_para= as.numeric(stringr::str_sub(tree_depth,-1,-1))
  
  #wide to long for MATRIX!
  dm <- reshape2::melt(d, id.vars = names(dimnames(d)), value.name = "value")
  
  dm$depth = rep(1:max.depth_para, p)
  
  #the line below is wrong code because "sort" function changed the variable order, the split frequency doesn't match variables anymore!!!
  #wrong code: dm <- data.frame(variable = sort(rep(names(d), nrow(d))), value = as.vector(as.matrix(d)), depth = rep(1:max.depth_para, p))
  
  for(i in 1:max.depth_para){
    tot.depth <- sum(dm[dm$depth == i,]$value)
    dm[dm$depth == i,]$value <- dm[dm$depth == i,]$value / tot.depth
  }
  

  g2 <- ggplot(dm, aes(x=variable, 
                       y=-depth, 
                       fill = value)) +
    geom_tile() +
    xlab("Variable") +
    ylab("Depth") +
    viridis::scale_fill_viridis(limits=c(0,1),discrete=FALSE) + 
    ggtitle( title_para,
             #paste(#tree_depth, CF, iteration_char  
                    #truth_description
                    #), 
             subtitle =   subtitle_para #paste(tree_depth, CF, iteration_char)
             ) + 
    theme(
      panel.background = element_blank(),
      plot.title    = element_text(color = "black", size=18, face="bold", # hjust = 0.5, 
                                   vjust=5, 
                                   margin =ggplot2:: margin(0.5, 0, 0, 0, "cm")),
      plot.subtitle = element_text(color = "black", size=16, face="bold" ),
      axis.title.x  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
      axis.title.y  = element_text(color = "black", size=14, face="plain",  hjust=0.5),
      plot.caption  = element_text(color = "black", size=12, face="plain"),
      axis.text.x   = element_text(color = "black", size=14, face="plain", angle = 0),
      axis.text.y   = element_text(color = "black", size=14, face="plain", angle = 0)) +
    scale_y_continuous(breaks=breaks_para, expand = c(0, 0), labels = label_para) +  
    scale_x_discrete(expand = c(0, 0)) +
    theme(legend.position = legend_para) +
    geom_text(aes(label = round(value, 2)), size = size_para) 

    return(g2)
  
}


#---------------------------------------------------------------
#Generating figure summary of HTE P-value for D4,D3,D2 iCF
#---------------------------------------------------------------

#' Function that produce P-value for miCF
#' @param df_D4 depth 4 CF 
#' @param df_D3 depth 3 CF
#' @param df_D2 depth 2 CF
#' @param type_i from all raw CF ("all") or those CF with best trees with the same structure as those majority voted best trees ("majority")
#' @param truth_describe describe the truth
#' @param iflegend if add legend
#' 
#' @return the dataframe of residual & MSE
#' 
#' @export
#' 
#Example: all_P_HTE   <- GG_HTE_P(HTE_P_iCF_D4,   HTE_P_iCF_D3,   HTE_P_iCF_D2,   "all",    truth_description, "nolegend")



#lab over

GG_HTE_P <- function(df_D4, df_D3, df_D2, type_i, truth_describe, iflegend){
  #stack P-value of different depth of iCF
  HTE_P_iCF_D432 <-rbind(df_D4, df_D3, df_D2)
  #make a factor to make sure y-aixs shown as D4 then D3 then D2
  HTE_P_iCF_D432$tree_depth_f <- factor(HTE_P_iCF_D432$tree_depth, levels =c("D4","D3","D2"), labels = c("D4 iCF","D3 iCF","D2 iCF"))
  #obtain depth of iCF with the maxium occurence (i.e. Stability) 
  df_count<- HTE_P_iCF_D432 %>% dplyr::count(tree_depth_f) %>% slice(which.max(n)) #maybe reload dplyr as it is affected by plyr
  HTE_P_iCF_D432$rawP <- HTE_P_cf.raw
  
  if (type_i == "majority") {
    title_para=paste0(#"iCF from voted best tree for ",
                      truth_describe)
    subTITLE = paste0("No. of forests for D4, D3, and D2 iCF: ", nrow(df_D4),", ", nrow(df_D3),", and ", nrow(df_D2), ", respectively" ) 
  } else  {
    title_para=paste0(#"iCF for ",
                      truth_describe)
    subTITLE = paste0("No. of forests: ", nrow(df_D4) ) 
    
  }
  
  scale_para = df_count$n/iterationNo
  #if (df_count$n/iterationNo < 0.15 ) {
  #  scale_para=0.1
  #} else  {
  #  scale_para=1
  #}
  
  
  if(iflegend=="nolegend"){
    legend_para = "none"}
  else {legend_para = "right"}
  
if(mean( round(  HTE_P_iCF_D432$rawP) )<0.001){
  Pvalue_para = "< 0.001"
  caption_para = "The raw CF P-value "
} else {
  Pvalue_para = mean( round(  HTE_P_iCF_D432$rawP ,3) )
  caption_para = "The raw CF P-value = "
  
}
  
  theme_set(theme_minimal())
  P_ridge <- ggplot(HTE_P_iCF_D432[ which(HTE_P_iCF_D432$HTE_P_cf >= 0), ],  
                    aes(x =HTE_P_cf, y=forcats::fct_rev(tree_depth_f), fill = stat(x)) ) +
    geom_density_ridges_gradient(scale =scale_para, #A scaling factor to scale the height of the ridgelines relative to the spacing between them. 
                                                    #A value of 1 indicates that the maximum point of any ridgeline touches the baseline right above, 
                                                    #assuming even spacing between baselines.
                                 size = 0.3, 
                                 rel_min_height = 0 #Lines with heights below this cutoff will be removed. 
                                 ) +
    scale_fill_viridis_c(name = "P-value", option = "D", begin = 0, end = 1,limits=c(0,1) ) +
    labs(title = title_para, 
         subtitle = , #subTITLE ,
         x="P-value", 
         y="Density" ,
         caption = paste0(caption_para, Pvalue_para, "(dotted line)") 
        ) +
    theme(
      plot.title    = element_text(color="black", size=18, face="bold"),
      plot.subtitle = element_text(color="black", size=16, face="bold"),
      axis.title.x  = element_text(color="black", size=14, face="plain"),
      axis.title.y  = element_text(color="black", size=14, face="plain"),
      plot.caption  = element_text(color="black", size=12, face="plain"),
      axis.text.x   = element_text(color="black", size=14, face="plain", angle=0),
      axis.text.y   = element_text(color="black", size=14, face="plain", angle=0))+
   #don't activate: geom_vline(aes(xintercept =   HTE_P_iCF_D432$rawP ),col='gray45', size=0.6 ,linetype="dashed") + #doesn't work when merging ggplot as the lastest vline will be applied to all other ggplots
    theme(legend.position = legend_para)
  return(P_ridge)
}


#-------------------------------
# Variable Importance
#-------------------------------
GG_VI <- function(varimp_df, title){
  VarImp <- varimp_df %>%
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


#' Function that shows subgroup Delta Y by subgrooup decisions obtained from each depth of forest (oneCF, iCF, or iCFv) 
#' @param model_df the df of coefficient (i.e. Delta Y) and its CI  
#' @methodL forest type: oneCF, iCF, or iCFv
#' @legend_para whether show legend of gradient fill, "nolgend" or "right"
#' @splitter splitter type
#' 
#' @return The ggplot
#' 
#' @export
#' 
GG_SG_DELTA_Y <- function(model_df, tree_depth, legend_para, splitter){
  
  if (model_df$Sub_def[1]=="W <100"){model_df$Sub_def="No Heterogeneous Subgroups"}
  
  if (splitter=="allBinary"){
  #first deal with the non-last subgroup definition
  model_df <- model_df %>% dplyr::mutate(Sub_def =   gsub ("<=0 &", "=0 &",Sub_def) )%>%
                           dplyr::mutate(Sub_def =   gsub (">0 &", "=1 &", Sub_def) )
  Sub_def_L <- as.list(model_df$Sub_def) 
  #str_sub doesn't work in column of dataframe or a list, thus have to use paste0 Fx to combine first and second part
  Sub_def_L_last_3 <- lapply(Sub_def_L, function(E) if(str_sub(E, -3, -1)=="<=0"){paste0( str_sub(E,1,-4), stringr::str_sub(E, -3,-1) <- "=0" ) } else{E} )
  
  Sub_def_L_last_3_2 <- lapply(Sub_def_L_last_3, function(E) if(str_sub(E, -2, -1)==">0"){paste0( str_sub(E,1,-3), stringr::str_sub(E, -2,-1) <- "=1" ) } else{E} )
  
  Sub_def_v = unlist(Sub_def_L_last_3_2, use.names=FALSE) 
  
  model_df <- model_df %>% mutate(Sub_def_v = Sub_def_v)
  
  }

          if(tree_depth=="D2" | tree_depth=="D3" | tree_depth=="D4"){
    title_para = ""
  }  else if (tree_depth=="truth") {
    title_para = paste(truth_description)
    
  } 

g3<- ggplot(model_df, aes(Estimate, column_label, color= )) + 
  geom_point(aes(x=Estimate, color = Estimate), size=5) + 
  geom_linerange(aes(xmin = conf.low, xmax=conf.high, color = Estimate),   size=1) + 
  geom_text(aes(label=Sub_def_v), vjust = 0, nudge_y = 0.2, size=2.5) +
  scale_color_viridis(option = "D", name = expression(paste(Delta, "Y")))+
  ggtitle(title_para) + 
  labs(subtitle= ifelse(tree_depth=="truth", paste("Truth"), paste(tree_depth, "iCF", sep = " ")), 
       caption = ifelse(tree_depth=="truth", paste(""),      paste0("Stability = ", mean(model_df$stability) ))  
       
         ) +
  ylab("Subgroup") + xlab(expression(paste("Subgroup ", Delta, "Y"))) +
  theme_bw() +
  theme(
    plot.title    = element_text(color = "black", size=16.5, face="bold", # hjust = 0.5, 
                                                                       vjust=5, 
                                                                       margin =ggplot2:: margin(0.5, 0, 0, 0, "cm")),
    plot.subtitle = element_text(color="black", size=16, face="bold"),
    axis.title.x  = element_text(color="black", size=14, face="plain"),
    axis.title.y  = element_text(color="black", size=14, face="plain"),
    plot.caption  = element_text(color="black", size=14, face="plain"),
    axis.text.x   = element_text(color="black", size=14, face="plain", angle=0),
    axis.text.y   = element_text(color="black", size=14, face="plain", angle=0))  + 
  scale_y_discrete(limits = unique(rev(model_df$column_label))) + 
  theme(legend.position = legend_para #"right"
        #"nolegend"
  )
return(g3)
}
