#################################################----------------------------------------------
#################################################----------------------------------------------
# II. TRUTH (TRUE TREE, TRUE SUBGROUPS) FOR INTERACTIONS
###################################CC4F419C|##############----------------------------------------------
#################################################----------------------------------------------
#' This function MAKE TRUE TREE DATAFRAME FOR DIFFERENT SCENARIOS,
#' "tree_true" is less important "tree_true_subgroup" as different tree structures may lead to the same subgroup decision
#' We make the "tree_true" dataframe so that "tree_true_subgroup" can be generated
#' bo4= binary:ordinary with 4 level
#' bo3= binary:ordinal with 3 level
#' 1o3= 1 ordinal 3 levels 
#' 1o4= 1 ordinal 4 levels 
#' @param intTRUE scenario name which provide information on n-way/antigonism or synergistic/moderate or strong/binary, categorical, or continuous, e.g."2W_Ant_Mo_1b"
#' 
#' @return The subsetted list that all are identical to reference list.
#' 
#' @export
#'
TRUTH <- function(intTRUE){
  if (	substr(intTRUE, 1, 2) == "5W" && stringr::str_sub(intTRUE, -2, -1) == "4b" ) {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07","node-08","node-09","node-10","node-11","node-12","node-13","node-14","node-15", "node-16","node-17","node-18","node-19","node-20","node-21","node-22","node-23","node-24","node-25","node-26","node-27","node-28","node-29","node-30","node-31"), 
                            "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE" , "FALSE"  , "FALSE" , "FALSE" , "FALSE" , "FALSE"  ,"FALSE"  ,"FALSE" , "FALSE" , "FALSE" , "FALSE" , "FALSE" ,"FALSE"  ,  "TRUE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   , "TRUE" ,"TRUE"    ,"TRUE"   ,"TRUE"   ,"TRUE"   , "TRUE"  , "TRUE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ),
                            "left_child"     = c("02"     ,"04"     ,"06"     ,"08"     ,"10"    ,  "12"    , "14"    , "16"     ,"18"     ,"20"     ,"22"     ,"24"    , "26"    , "28"    , "30"    ,  "NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"    ,"NA"      ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                            "right_child"    = c("03"     ,"05"     ,"07"     ,"09"     ,"11"    ,  "13"    , "15"    , "17"     ,"19"     ,"21"     ,"23"     ,"25"    , "27"    , "29"    , "31"    ,  "NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"    ,"NA"      ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                            "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"X8"     ,"X8"    ,  "X8"    , "X8"    , "X9"     ,"X9"     ,"X9"     ,"X9"     ,"X9"    , "X9"    , "X9"    , "X9"    ,  "NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"    ,"NA"      ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                            "split_value"    = c("0"      ,"0"      ,"0"      ,"0"      ,"0"     ,  "0"     , "0"     , "0"      ,"0"      ,"0"      ,"0"      ,"0"     , "0"     , "0"     ,  "0"    ,  "NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"    ,"NA"      ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                            stringsAsFactors=FALSE)
    
    tree_depth_true  <<-"D5"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- tree_true[ 1:3, ]
    tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
    truth_description <<- "5-way INT of W, X3, X1, X8, X9"
    truth_INT <<- c("W:X1:X3:X8:X9", 
                 "W:X3:X8:X9", "W:X1:X3:X9", "W:X1:X3:X8", 
                 "W:X1:X3", "W:X1:X8","W:X1:X9", "W:X3:X8",  "W:X3:X9", "W:X8:X9",
                 "W:X3", "W:X1", "W:X8", "W:X9")
    vars_catover2 <<- NA
  } else if (	substr(intTRUE, 1, 2) == "4W" && stringr::str_sub(intTRUE, -2, -1) == "3b" ) {
  tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07","node-08","node-09","node-10","node-11","node-12","node-13","node-14","node-15"), 
                          "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE" ,"FALSE" ,"FALSE" ,"FALSE" ,"FALSE" ,"TRUE"  ,"TRUE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ),
                          "left_child"     = c("02"     ,"04"     ,"06"     ,"08"     ,"10"    ,"12"    ,"14"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "right_child"    = c("03"     ,"05"     ,"07"     ,"09"     ,"11"    ,"13"    ,"15"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"X8"     ,"X8"    ,"X8"    ,"X8"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "split_value"    = c("0"      ,"0"      ,"0"      ,"0"      ,"0"     ,"0"     ,"0"     ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          stringsAsFactors=FALSE)
  
  tree_depth_true  <<-"D4"
  tree_true_N1     <<- tree_true[ 1, ]
  tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
  tree_true_N123   <<- tree_true[ 1:3, ]
  tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
  truth_description <<- "4-way INT of W, X3, X1, X8"
  truth_INT         <<-c("W:X1:X3:X8", "W:X1:X3", "W:X1:X8", "W:X3:X8", "W:X3", "W:X1", "W:X8")
  vars_catover2     <<- NA
} else if (	substr(intTRUE, 1, 2) == "4W" && stringr::str_sub(intTRUE, -4, -1) == "bbo4" ) {
  tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07","node-08","node-09","node-10","node-11","node-12","node-13","node-14","node-15"), 
                          "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE" ,"FALSE" ,"FALSE" ,"FALSE" ,"FALSE" ,"TRUE"  ,"TRUE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ),
                          "left_child"     = c("02"     ,"04"     ,"06"     ,"08"     ,"10"    ,"12"    ,"14"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "right_child"    = c("03"     ,"05"     ,"07"     ,"09"     ,"11"    ,"13"    ,"15"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"X2"     ,"X2"    ,"X2"    ,"X2"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "split_value"    = c("0"      ,"0"      ,"0"      ,"1"      ,"1"     ,"1"     ,"1"     ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          stringsAsFactors=FALSE)
  
  tree_depth_true  <<-"D4"
  tree_true_N1     <<- tree_true[ 1, ]
  tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
  tree_true_N123   <<- tree_true[ 1:3, ]
  tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
  truth_description <<- "4-way INT of W, X3, X1, X2(4l)"
  truth_INT        <<- c("W:X1:X2:X3", "W:X1:X3", "W:X1:X2", "W:X2:X3", "W:X3", "W:X1", "W:X2")
  vars_catover2    <<- "X2"
} else if (	substr(intTRUE, 1, 2) == "4W" && stringr::str_sub(intTRUE, -3, -1) == "bbl" ) {
  tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07","node-08","node-09","node-10","node-11","node-12","node-13","node-14","node-15"), 
                          "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE" ,"FALSE" ,"FALSE" ,"FALSE" ,"FALSE" ,"TRUE"  ,"TRUE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ),
                          "left_child"     = c("02"     ,"04"     ,"06"     ,"08"     ,"10"    ,"12"    ,"14"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "right_child"    = c("03"     ,"05"     ,"07"     ,"09"     ,"11"    ,"13"    ,"15"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"X2"     ,"X2"    ,"X2"    ,"X2"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          "split_value"    = c("0"      ,"0"      ,"0"      ,"0.5"      ,"0.5"     ,"0.5"     ,"0.5"     ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                          stringsAsFactors=FALSE)
  
  tree_depth_true  <<-"D4"
  tree_true_N1     <<- tree_true[ 1, ]
  tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
  tree_true_N123   <<- tree_true[ 1:3, ]
  tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
  truth_description <<- "4-way INT of W, X3, X1, X2(l)"
  truth_INT         <<-c("W:X1:X2:X3", "W:X1:X3", "W:X1:X2", "W:X2:X3", "W:X3", "W:X1", "W:X2")
  vars_catover2     <<- NA
}else if (  substr(intTRUE, 1, 2) == "3W" && stringr::str_sub(intTRUE, -2, -1) == "2b") {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07"), 
                            "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"  ),
                            "left_child"     = c("02"     ,"04"     ,"06"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "right_child"    = c("03"     ,"05"     ,"07"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "split_value"    = c("0"      ,"0"      ,"0"      ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            stringsAsFactors=FALSE)
    tree_depth_true<<-"D3"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- tree_true[ 1:3, ]
    tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
    truth_description <<- "3-way INT of W, X3, X1"
    truth_INT <<-c("W:X1:X3", "W:X3", "W:X1")
    vars_catover2 <<- NA
  } else if (	substr(intTRUE, 1, 2) == "3W" && stringr::str_sub(intTRUE, -3, -1) == "bo4" ) {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07"), 
                            "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE"  ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  ),
                            "left_child"     = c("02"     ,"04"     ,"06"     ,"NA"     ,"NA"     ,"NA"    ,"NA"    ),
                            "right_child"    = c("03"     ,"05"     ,"07"     ,"NA"     ,"NA"     ,"NA"    ,"NA"    ),
                            "split_variable" = c("X3"     ,"X2"     ,"X2"     ,"NA"     ,"NA"     ,"NA"    ,"NA"    ),
                            "split_value"    = c("0"      ,"1"      ,"1"      ,"NA"     ,"NA"     ,"NA"    ,"NA"    ),
                            stringsAsFactors=FALSE)
    tree_depth_true<<-"D3"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- tree_true[ 1:3, ]
    tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
    truth_description <<- "3-way INT of W, X3, X2(4l)"
    truth_INT         <<-c("W:X2:X3", "W:X3", "W:X2")
    vars_catover2     <<- "X2"
  } else if (  substr(intTRUE, 1, 2) == "3W" && stringr::str_sub(intTRUE, -3, -1) == "bo3") {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07"), 
                            "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"  ),
                            "left_child"     = c("02"     ,"04"     ,"06"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "right_child"    = c("03"     ,"05"     ,"07"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "split_variable" = c("X3"     ,"X2"     ,"X2"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "split_value"    = c("0"      ,"1"      ,"1"      ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            stringsAsFactors=FALSE)
    tree_depth_true  <<-"D3"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- tree_true[ 1:3, ]
    tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
    truth_description <<- "3-way INT of W, X3, X2(3l)"
    truth_INT         <<-c("W:X2:X3", "W:X3", "W:X2")
    vars_catover2     <<- "X2"
  } else if (  substr(intTRUE, 1, 2) == "3W" && stringr::str_sub(intTRUE, -2, -1) == "bl") {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07"), 
                            "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"  ),
                            "left_child"     = c("02"     ,"04"     ,"06"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "right_child"    = c("03"     ,"05"     ,"07"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "split_variable" = c("X3"     ,"X2"     ,"X2"     ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            "split_value"    = c("0"      ,"0.5"      ,"0.5"      ,"NA"     ,"NA"     ,"NA"     ,"NA"    ),
                            stringsAsFactors=FALSE)
    tree_depth_true  <<-"D3"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- tree_true[ 1:3, ]
    tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
    truth_description <<- "3-way INT of W, X3, X2(l)"
    truth_INT         <<-c("W:X2:X3", "W:X3", "W:X2")
    vars_catover2     <<- NA
  } else if (  substr(intTRUE, 1, 2) == "2W" && stringr::str_sub(intTRUE, -2, -1) == "1b" ) {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03"), 
                            "is_leaf"        = c("FALSE"  ,"TRUE"   ,"TRUE"  ),
                            "left_child"     = c("02"     ,"NA"     ,"NA"    ),
                            "right_child"    = c("03"     ,"NA"     ,"NA"    ),
                            "split_variable" = c("X3"     ,"NA"     ,"NA"    ),
                            "split_value"    = c("0"      ,"NA"     ,"NA"    ),
                            stringsAsFactors=FALSE)   
    tree_depth_true  <<-"D2"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- "NA"
    tree_true_N123_r <<- "NA"
    truth_description <<- "2-way INT of W, X3"
    truth_INT        <<- c("W:X3") 
    vars_catover2     <<- NA
  } else if (  substr(intTRUE, 1, 2) == "2W" && stringr::str_sub(intTRUE, -2, -1) == "1l" ) {
      tree_true <<- data.frame("node"           = c("node-01","node-02","node-03"), 
                              "is_leaf"        = c("FALSE"  ,"TRUE"   ,"TRUE"  ),
                              "left_child"     = c("02"     ,"NA"     ,"NA"    ),
                              "right_child"    = c("03"     ,"NA"     ,"NA"    ),
                              "split_variable" = c("X2"     ,"NA"     ,"NA"    ),
                              "split_value"    = c("0.5"      ,"NA"     ,"NA"    ),
                              stringsAsFactors=FALSE)   
      tree_depth_true  <<-"D2"
      tree_true_N1     <<- tree_true[ 1, ]
      tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
      tree_true_N123   <<- "NA"
      tree_true_N123_r <<- "NA"
      truth_description <<- "2-way INT of W, X2(l)"
      truth_INT         <<-c("W:X2")
      vars_catover2      <<- NA
  } else if (  substr(intTRUE, 1, 2) == "2W" && stringr::str_sub(intTRUE, -3, -1) == "1o3" ) {
      tree_true <<- data.frame("node"           = c("node-01","node-02","node-03"), 
                              "is_leaf"        = c("FALSE"  ,"TRUE"   ,"TRUE"  ),
                              "left_child"     = c("02"     ,"NA"     ,"NA"    ),
                              "right_child"    = c("03"     ,"NA"     ,"NA"    ),
                              "split_variable" = c("X2"     ,"NA"     ,"NA"    ),
                              "split_value"    = c("1"      ,"NA"     ,"NA"    ),
                              stringsAsFactors=FALSE)   
      tree_depth_true  <<-"D2"
      tree_true_N1     <<- tree_true[ 1, ]
      tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
      tree_true_N123   <<- "NA"
      tree_true_N123_r <<- "NA"
      truth_description <<- "2-way INT of W, X2(3l)"
      truth_INT         <<-c("W:X2")
      vars_catover2     <<- "X2"
  } else if (  substr(intTRUE, 1, 2) == "2W" && stringr::str_sub(intTRUE, -3, -1) == "1o4" ) {
    tree_true <<- data.frame("node"           = c("node-01","node-02","node-03"), 
                            "is_leaf"        = c("FALSE"  ,"TRUE"   ,"TRUE"  ),
                            "left_child"     = c("02"     ,"NA"     ,"NA"    ),
                            "right_child"    = c("03"     ,"NA"     ,"NA"    ),
                            "split_variable" = c("X2"     ,"NA"     ,"NA"    ),
                            "split_value"    = c("1"      ,"NA"     ,"NA"    ),
                            stringsAsFactors=FALSE)   
    tree_depth_true  <<-"D2"
    tree_true_N1     <<- tree_true[ 1, ]
    tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
    tree_true_N123   <<- "NA"
    tree_true_N123_r <<- "NA"
    truth_description <<- "2-way INT of W, X2(4l)"
    truth_INT         <<-c("W:X2")
    vars_catover2     <<- "X2"
  } else if (  substr(intTRUE, 1, 2) == "2W" && stringr::str_sub(intTRUE, -2, -1) == "2b" ) {
      tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07"), 
                              "is_leaf"        = c("FALSE" ,"FALSE"   ,"FALSE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ),
                              "left_child"     = c("02"     ,"04"     ,"06"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              "right_child"    = c("03"     ,"05"     ,"07"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              "split_value"    = c("0"      ,"0"      ,"0"      ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              stringsAsFactors=FALSE)
      tree_depth_true  <<-"D3"
      tree_true_N1     <<- tree_true[ 1, ]
      tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
      tree_true_N123   <<- tree_true[ 1:3, ]
      tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
      
      truth_description <<- "Two 2-way INTs of W, X3, X1"
      truth_INT         <<- c("W:X3", "W:X1")
      vars_catover2     <<- NA
      
  } else if (  substr(intTRUE, 1, 2) == "2W" && stringr::str_sub(intTRUE, -2, -1) == "3b" ) {
      tree_true <<- data.frame("node"           = c("node-01","node-02","node-03","node-04","node-05","node-06","node-07","node-08","node-09","node-10","node-11","node-12","node-13","node-14","node-15"), 
                              "is_leaf"        = c("FALSE"  ,"FALSE"  ,"FALSE"  ,"FALSE"  ,"FALSE"  ,"FALSE"  ,"FALSE" ,"TRUE"  ,"TRUE"  ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ,"TRUE"   ),
                              "left_child"     = c("02"     ,"04"     ,"06"     ,"08"     ,"10"     ,"12"     ,"14"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              "right_child"    = c("03"     ,"05"     ,"07"     ,"09"     ,"11"     ,"13"     ,"15"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              "split_variable" = c("X3"     ,"X1"     ,"X1"     ,"X8"     ,"X8"     ,"X8"     ,"X8"    ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              "split_value"    = c("0"      ,"0"      ,"0"      ,"0"      ,"0"      ,"0"      ,"0"     ,"NA"    ,"NA"    ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ,"NA"     ),
                              stringsAsFactors=FALSE)
      tree_depth_true  <<-"D4"
      tree_true_N1     <<- tree_true[ 1, ]
      tree_true_N1_r   <<- tree_true[ 1, 1:5 ]
      tree_true_N123   <<- tree_true[ 1:3, ]
      tree_true_N123_r <<- tree_true[ 1:3, 1:5 ]
      truth_description <<- "Three 2-way INTs of W, X3, X1, X8" #W-X3, W-X1, W-X8
      truth_INT         <<- c("W:X3", "W:X1", "W:X8")
      vars_catover2     <<- NA
  } else if (  substr(intTRUE, 1, 2) == "0W" ) {
        tree_true        <<- "NA" 
        tree_depth_true  <<- "NA" 
        tree_true_N1     <<- "NA" 
        tree_true_N1_r   <<- "NA" 
        tree_true_N123   <<- "NA" 
        tree_true_N123_r <<- "NA" 
        truth_description <<- "No HTE"
        truth_INT         <<- "NA"  
        vars_catover2 <<- NA
  } else if (   substr(intTRUE, 1, 2) == "Un") {
    tree_true         <<- "Unknown" 
    tree_depth_true   <<- "Unknown"  
    tree_true_N1      <<- "Unknown"  
    tree_true_N1_r    <<- "Unknown"  
    tree_true_N123    <<- "Unknown"  
    tree_true_N123_r  <<- "Unknown"  
    truth_description <<- "Unknown" 
    truth_INT         <<- "Unknown"  
    #vars_catover2     <<- "Assign vars_catover2 Manually"
  }


#relaxed true tree, i.e. ignoring split values
if (  substr(intTRUE, 1, 2) == "0W" ) {
  tree_true_r        <<- "NA" 
  tree_true_subgroup <<- "NA"  #ohter scearios calculated below
  truth_INT          <<- "NA"  
  
} else if (substr(intTRUE, 1, 2) != "Un") {
  #remove "split_value" column from tree_true:
  tree_true_r <<- tree_true[c(-6)] 
  #----------------------------------------------------------------------------------------------
  #obtain subgroup decision for the truth, taking advantage of PRE_MAJORITY_SUBGROUP function
  #----------------------------------------------------------------------------------------------
  #--------- the following 2 lines are intermediate steps--------------------------
  #make a list first so that PRE_MAJORITY_SUBGROUP can run properly, keep this list format.
  tree_true_subgroup_0   <<- PRE_MAJORITY_SUBGROUP ( list(tree_true)  ) 
  #drop other covariates
  tree_true_subgroup_key <<-lapply(tree_true_subgroup_0, function(df) subset(df, select=c("subgroupID", "subgroup")))
  #--------------------------------------------------------------------------------
  tree_true_subgroup <<- do.call(rbind.data.frame, tree_true_subgroup_key) #covert a list to a dataframe
  
} else if ( substr(intTRUE, 1, 2) == "Un" ){
  
  tree_true_r            <<-  "Unknown"  
  tree_true_subgroup_0   <<- "Unknown"   
  tree_true_subgroup_key <<- "Unknown" 
  tree_true_subgroup     <<- "Unknown" 
}
truth.list <- list(tree_true          = tree_true,
                   tree_depth_true    = tree_depth_true,
                   tree_true_N1       = tree_true_N1,
                   tree_true_N1_r     = tree_true_N1_r,
                   tree_true_N123     = tree_true_N123,
                   tree_true_N123_r   = tree_true_N123_r,
                   tree_true_r        = tree_true_r,
                   tree_true_subgroup = tree_true_subgroup,
                   truth_description  = truth_description,
                   truth_INT          = stringr::str_sort( truth_INT)
                   #vars_catover2      = vars_catover2 
                   )
return(truth.list)
                   
}






#' function to assess performance of different tree-based subgrouping methods 
#' @param Deci_SG subgroup decision 
#' @param iterationNo distinguish between oneCF and iCF
#'  
#' @return final subgroup decision G_iCF
#' 
#' @export
#' 
#'

TREE_PERFORMANCE <- function(Deci_SG, iterationNo){
  #performance 
  DltY_iCF_te      = DltY_DATA(Deci_SG,       Test_ID, "predict" )
  MSE_iCF_te       <- MSE_DltY(DltY_true_te, DltY_iCF_te)
  MSE_ate_iCF_te   <- as.numeric(MSE_iCF_te$MSE_ate)
  MSE_att_iCF_te   <- as.numeric(MSE_iCF_te$MSE_att)
  
  #DISCOVERYRATE function take care of the "NA" decision
  if(iterationNo >1){
    disco_SG_iCF       = DISCOVERYRATE(Deci_SG ,                 truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")$value
    disco_INT_iCF      = DISCOVERYRATE(SG2INT(Deci_SG),          truth_description, tree_true_subgroup, truth_INT, "iCF", "interaction")$value
    #disco_SG_iCF.act.aa    = DISCOVERYRATE(Deci_SG,                  truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")$sg_Nacc_Nall
    #disco_SG_iCF.act.aat   = DISCOVERYRATE(Deci_SG,                  truth_description, tree_true_subgroup, truth_INT, "iCF", "subgroup")$sg_Nacc_NallT
    disco_CATEmax_iCF      = DISCOVERYRATE(CATE_MAX_SG(DltY_iCF_te), truth_description, tree_true_subgroup, CATE_max_true_te, "iCF", "CATEmax")$value
    
  } else if (iterationNo==1) {
    disco_SG_iCF       = DISCOVERYRATE(Deci_SG ,                  truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$value
    disco_INT_iCF      = DISCOVERYRATE(SG2INT(Deci_SG),           truth_description, tree_true_subgroup, truth_INT, "oneCFb", "interaction")$value
    #disco_SG_iCF.act.aa    = DISCOVERYRATE(Deci_SG,                   truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$sg_Nacc_Nall
    #disco_SG_iCF.act.aat   = DISCOVERYRATE(Deci_SG,                   truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$sg_Nacc_NallT
    disco_CATEmax_iCF      = DISCOVERYRATE(CATE_MAX_SG(DltY_iCF_te),  truth_description, tree_true_subgroup, CATE_max_true_te, "oneCFb", "CATEmax")$value
    
  }  
  return(list(disco_SG_iCF   = disco_SG_iCF ,
              disco_INT_iCF  = disco_INT_iCF, 
              #disco_SG_AIC_iCF.act.aa = disco_SG_AIC_iCF.act.aa,
              #disco_SG_AIC_iCF.act.aat  = disco_SG_AIC_iCF.act.aat,
              #DltY_iCF_te       = DltY_iCF_te,
              disco_CATEmax_iCF = disco_CATEmax_iCF,
              MSE_ate_iCF_te    = MSE_ate_iCF_te,
              MSE_att_iCF_te    = MSE_att_iCF_te)
  )  
  
}
