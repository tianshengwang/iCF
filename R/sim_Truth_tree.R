#################################################----------------------------------------------
#################################################----------------------------------------------
# II. TRUTH (TRUE TREE, TRUE SUBGROUPS) FOR INTERACTIONS
###################################CC4F419C|##############----------------------------------------------
#################################################----------------------------------------------
#' TRUTH
#' 
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
    disco_CATEmax_iCF      = DISCOVERYRATE(CATE_MAX_SG(DltY_iCF_te), truth_description, tree_true_subgroup, CATE_max_true_te, "iCF", "CATEmax")$value
    
  } else if (iterationNo==1) {
    disco_SG_iCF       = DISCOVERYRATE(Deci_SG ,                  truth_description, tree_true_subgroup, truth_INT, "oneCFb", "subgroup")$value
    disco_INT_iCF      = DISCOVERYRATE(SG2INT(Deci_SG),           truth_description, tree_true_subgroup, truth_INT, "oneCFb", "interaction")$value
    disco_CATEmax_iCF      = DISCOVERYRATE(CATE_MAX_SG(DltY_iCF_te),  truth_description, tree_true_subgroup, CATE_max_true_te, "oneCFb", "CATEmax")$value
    
  }  
  return(list(disco_SG_iCF   = disco_SG_iCF ,
              disco_INT_iCF  = disco_INT_iCF, 
              disco_CATEmax_iCF = disco_CATEmax_iCF,
              MSE_ate_iCF_te    = MSE_ate_iCF_te,
              MSE_att_iCF_te    = MSE_att_iCF_te)
  )  
  
}



#'CONVERT_tree_IT2CF
#'
#' This function is convert "interaction tree" into the tree format of a causal forest
#' @param btree_it interaction tree
#' 
#' @return  the casual forest format for interaction tree 
#'
#' @export
#' 

CONVERT_tree_IT2CF <- function(btree_it) {
  btree_cf <- btree_it  %>%  
    mutate (is_leaf        = ifelse(is.na(vname)==TRUE , TRUE , FALSE))  %>%
    mutate (left_child     = ifelse(is_leaf=="FALSE", paste0(node,1), "NA")) %>%
    mutate (right_child    = ifelse(is_leaf=="FALSE", paste0(node,2), "NA")) %>%
    mutate (split_variable = ifelse(is.na(vname)==FALSE, paste0(as.character(vname)), "NA" ) ) %>%
    mutate (split_value    = ifelse(is.na(best.cut)==FALSE, round(as.numeric(as.character(best.cut)),1), "NA") ) %>%
    mutate (samples        = ifelse(is_leaf=="TRUE", n, "NA") ) %>%
    mutate (avg_Y          = 0) %>%
    mutate (avg_W          = 0) %>%
    mutate (b              = 1) %>%
    mutate (HTE_P_cf       = 0) %>%
    subset (select=  c(node, is_leaf, left_child, right_child, split_variable, split_value, samples, avg_Y, avg_W, b, HTE_P_cf)) %>%
    mutate (node = paste("node", node, sep = '-') ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -2, -2) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -3, -3) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -4, -4) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -5, -5) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (left_child  = ifelse(nchar(left_child)  == 2 & left_child  !="NA", paste0("000", left_child),  left_child ) ) %>%
    mutate (left_child  = ifelse(nchar(left_child)  == 3 & left_child  !="NA", paste0("00" , left_child),  left_child ) ) %>%
    mutate (left_child  = ifelse(nchar(left_child)  == 4 & left_child  !="NA", paste0("0" ,  left_child),  left_child ) ) %>%
    mutate (right_child = ifelse(nchar(right_child) == 2 & right_child !="NA", paste0("000", right_child), right_child ) ) %>%
    mutate (right_child = ifelse(nchar(right_child) == 3 & right_child !="NA", paste0("00" , right_child), right_child ) ) %>%
    mutate (right_child = ifelse(nchar(right_child) == 4 & right_child !="NA", paste0("0" ,  right_child), right_child ) ) %>%
    arrange(node) 
  return( as.data.frame(btree_cf) )       
}


#'CONVERT_tree_stima2CF
#'
#' This function is convert "STIMA tree" into the tree format of a causal forest
#' @param btree_stima STIMA tree
#' 
#' @return the casual forest format for STIMA tree 
#'
#' @export
#'
CONVERT_tree_stima2CF <- function(btree_stima) {
  btree_cf <- btree_stima %>% 
    mutate (node     = 1:nrow(btree_stima)) %>%
    mutate (node     = paste("node",node, sep = "-" )) %>%
    mutate (node = ifelse( stringr::str_sub(node,  -2, -2) =="-", sub("-([0-9])", "-0\\1", node), node ) ) %>%
    mutate (is_leaf     = ifelse(nchar(Region)!=0 , TRUE , FALSE)) %>%
    mutate ( value_1    =  lead(Splitpoint,n=1)   )%>%  
    mutate ( value_2    =  lead(Splitpoint,n=2)   )%>% 
    mutate ( value_3    =  lead(Splitpoint,n=3)   )%>%  
    mutate ( value_4    =  lead(Splitpoint,n=4)   )%>% 
    mutate ( var_1      =  lead(Predictor,n=1)   )   %>%  
    mutate ( var_2      =  lead(Predictor,n=2)   )   %>%  
    mutate ( var_3      =  lead(Predictor,n=3)   )   %>%  
    mutate ( var_4      =  lead(Predictor,n=4)   )   %>% 
    mutate ( Sign_1          =  lead(Sign,n=1)   )   %>%  
    mutate ( Sign_2          =  lead(Sign,n=2)   )   %>%  
    mutate ( Sign_3          =  lead(Sign,n=3)   )   %>%  
    mutate ( Sign_4          =  lead(Sign,n=4)   )   %>% 
    mutate ( node_1     =  lead(node,n=1)   )   %>%  
    mutate ( node_2     =  lead(node,n=2)   )   %>%  
    mutate ( node_3     =  lead(node,n=3)   )   %>%  
    mutate ( node_4     =  lead(node,n=4)   )   %>% 
    mutate ( Region_1     =  lead(Region,n=1)   )   %>%  
    mutate ( Region_2     =  lead(Region,n=2)   )   %>%  
    mutate ( Region_3     =  lead(Region,n=3)   )   %>%  
    mutate ( Region_4     =  lead(Region,n=4)   )   %>% 
    mutate (left_child =
              ifelse(is_leaf=="TRUE", "NA",          
                     ifelse((Sign_1 == "<=" & value_1== value_2 & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, stringr::str_sub(node_1,-2,-1) , 
                            ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                    & value_2== value_3 & var_2==var_3 &  ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, stringr::str_sub(node_2,-2,-1),
                                   ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                           & value_1== value_2 & #value_3== value_4 & 
                                             var_1==var_2 #& var_3==var_4 
                                           & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, stringr::str_sub(node_3,-2,-1), NA
                                   ))))) %>%
    mutate (right_child =
              ifelse(is_leaf=="TRUE", "NA",
                     ifelse((Sign_1 == "<=" & value_1== value_2 & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, stringr::str_sub(node_2,-2,-1) , 
                            ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                    & value_2== value_3 & var_2==var_3 & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) ) ==TRUE, stringr::str_sub(node_3,-2,-1),
                                   ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                           & value_1== value_2 & value_3== value_4 & var_1==var_2 & var_3==var_4 & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, stringr::str_sub(node_4,-2,-1), NA
                                   ))))) %>%
    mutate (split_variable = 
              ifelse(is_leaf=="TRUE", "NA",          
                     ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"  
                             & value_1== value_2 #& value_3== value_4 
                             &( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, var_1, 
                            ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                    #& value_2== value_3 & var_2==var_3 
                                    & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) ) ==TRUE,    var_2,
                                   ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                           & value_1== value_2 & value_3== value_4 & var_1==var_2 & var_3==var_4 & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, var_3, NA
                                   ))))) %>%
    mutate (split_value = 
              ifelse(is_leaf=="TRUE", "NA",          
                     ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"  
                             & value_1== value_2 #& value_3== value_4 
                             & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) )==TRUE, round(as.numeric(as.character(value_1)),1), 
                            ifelse((Sign_1 == ">"  & Sign_2 == "<=" #& Sign_3 == ">"  & Sign_4 == "<="  
                                    #& value_2== value_3 & var_2==var_3 
                                    & ( nchar(Region_1)==0|is.na(Region_1) |nchar(Region_2)==0|is.na(Region_2)|nchar(Region_3)==0|is.na(Region_3)|nchar(Region_4)==0|is.na(Region_4) ) ) ==TRUE,   round(as.numeric(as.character(value_2)),1),
                                   ifelse((Sign_1 == "<=" & Sign_2 == ">"  #& Sign_3 == "<=" & Sign_4 == ">"   
                                           & value_1== value_2 & value_3== value_4 & var_1==var_2 & var_3==var_4 & (nchar(Region_1)>=2&nchar(Region_2)>=2&nchar(Region_3)>=2&nchar(Region_4)>=2) ) ==TRUE, round(as.numeric(as.character(value_3)),1), NA
                                   )))))  %>%
    
    mutate (samples        = ifelse(is_leaf=="TRUE", n, "NA") ) %>%
    mutate (avg_Y          = 0) %>%
    mutate (avg_W          = 0) %>%
    mutate (b              = 1) %>%
    mutate (HTE_P_cf       = 0) %>%
    subset (select=  c(node, is_leaf, left_child, right_child, split_variable, split_value, samples, avg_Y, avg_W, b, HTE_P_cf))
  
  return( as.data.frame(btree_cf) )       
}


#

#' N_SUBGROUP
#' Function check number of subgroups obtained, and is used in the DISCOVERYRATE function below
#' @param deci list
#' @param method_L vector of string
#' 
#' @return T of F
#' 
#' @export
#' 
N_SUBGROUP <- function(deci, method_L) {
  if (is.null (nrow(deci))==TRUE  & method_L == "tree_sg" ) {
    n_subgroup = 0
  }  else if (is.null (nrow(deci))==FALSE & method_L == "tree_sg") {
    n_subgroup = nrow(deci)
  }  else if (deci[1] =="NA"              & method_L == "non-tree_int" ){
    n_subgroup = 0
  }  else if (deci[1] !="NA"              & method_L == "non-tree_int" ){
    n_subgroup = length(deci)
  }
  return(n_subgroup)
}


#' DETECT_STRING
#'
#' Function check if vector string in 2nd position nested in the list in 1st position, will be used in DISCOVERYRATE function
#' https://community.rstudio.com/t/detect-presence-absence-of-elements-in-nested-lists/31515
#' @param your_list list
#' @param vector_strings vector of string
#' @return T of F
#' 
#' @export

DETECT_STRING <- function(your_list, vector_strings){
  lapply(your_list, function(x) {
    if(TRUE %in% str_detect(x, paste(vector_strings, collapse = "|"))){
      TRUE
    } else {FALSE}
  })
}

#' DISCOVERYRATE
#' 
#' Function that assess accuracy, false positive, and false negative,
#' which is used for assessing performance of difference SG/INT identification methods in simulation study of the iCF paper
#' Accuracy definition: 
#' all subgup decisions (for tree-based methods that explicitly give subgroup) or all interactions (for non-tree based methds) are accurately identified  
#' @param decision the final decision could be subgroup dataframe or interction vector
#' @param truth_description the description of truth
#' @param tree_true_subgroup the true subgroup
#' @param truth_INT the true interaction, when the metric==CATEmax, use true CATEmax subgroup instead of true interaction!!!
#' @param method different methods (tree-based or non-tree-based)
#' @param metric assesing perforamnce of methods by subgroup or interaction identification
#' 
#' @return value=score of accuracy assessment; sg_Nacc_Nall=% of true subgroups among identified subgroups; sg_Nacc_NallT=% of true subgroups among total true subgroups 
#' 
#' @export

DISCOVERYRATE <- function (decision, truth_description, tree_true_subgroup, truth_INT, method, metric){
  
  #if assessing performance by subgroup identification for tree-based method:
  if ( method %in% c("oneCFv", "oneCFb", "iCF", "iCFv", "IT", "VT")==T & metric=="subgroup"){
    if (is.null( dim(decision))== F ) {
      Deci_vec    <- dplyr::pull(decision,           subgroup)#subgroup here is a column's name
    } else {
      Deci_vec = "NA"
    }
    
    if (is.null( dim(tree_true_subgroup)) == F) {
      Truth_vec  <- dplyr::pull(tree_true_subgroup, subgroup)
    } else {
      Truth_vec = "NA"
    }
    #get N to judge false+ or false-
    N_true   <- N_SUBGROUP (tree_true_subgroup, "tree_sg")
    N_deci   <- N_SUBGROUP (decision,           "tree_sg")
  } else if ( metric=="interaction" | metric=="CATEmax") {
    #if assessing performance by interaction identification for ALL method:
    
    
    Deci_vec       <- decision
    Truth_vec      <- truth_INT
    N_true  <- N_SUBGROUP (truth_INT,      "non-tree_int")
    N_deci  <- N_SUBGROUP (decision,       "non-tree_int")
  }
  
  Deci_vec_L  <- as.list(Deci_vec)  #the same for both TREE and non-TREE method 
  Truth_vec_L <- as.list(Truth_vec) #the same for both TREE and non-TREE method
  
  value = ifelse (
    ###########################--------------------------       
    # I.1 ACCURATE
    ###########################--------------------------
    (
      #----------------------------------------------------                                       
      #I.1.1. identical for tree-based method
      #----------------------------------------------------  
      ( metric=="subgroup"                           & identical(decision, tree_true_subgroup)  ) |
        ((metric=="CATEmax" |  metric=="interaction" ) & identical(decision, truth_INT #when metric=="CATEmax", truth_INT parameter is true subgroup with max CATE rather than overall true interaction
        )  )  |
        #----------------------------------------------------                                       
      #I.1.2. 4-way INT
      #----------------------------------------------------
      (metric=="interaction" & truth_description == "4-way INT of W, X3, X1, X8" & 
         #I.1.2.1. the same ignoring order
         ( setequal(Deci_vec, truth_INT) == T  | 
             #I.1.2.2. as long as W:X1:X3:X8 is identified:
             ( "W:X1:X3:X8" %in% decision  & 
                 #all inconsistent interactions are one of the following:
                 all( decision[is.na(match(decision, truth_INT))] %in% c("W:X1:X3", "W:X1:X8", "W:X3:X8", "W:X1", "W:X3", "W:X8") )  
             )
         )
      )
      |
        (metric=="interaction" & (truth_description == "4-way INT of W, X3, X1, X2(4l)" |  #ordinal splitter 4 levels
                                    truth_description == "4-way INT of W, X3, X1, X2(3l)" |  #ordinal splitter 3 levels
                                    truth_description == "4-way INT of W, X3, X1, X2(l)") &  #cotinuous splitter
           #I.1.3.1. the same ignoring order
           ( setequal(Deci_vec, truth_INT) == T | 
               #I.1.3.2. as long as W:X1:X3 is identified:
               ("W:X1:X2:X3" %in% decision  & 
                  #all inconsistent interactions are one of the following:
                  all( decision[is.na(match(decision, truth_INT))] %in% c("W:X1:X3", "W:X1:X2", "W:X2:X3", "W:X1", "W:X3", "W:X2") )  
               ) 
           )
        )
      |
        #----------------------------------------------------                                       
      #I.1.3. 3-way INT
      #----------------------------------------------------
      (metric=="interaction" & truth_description == "3-way INT of W, X3, X1"     & #binary splitter
         #I.1.3.1. the same ignoring order
         ( setequal(Deci_vec, truth_INT) == T | 
             #I.1.3.2. as long as W:X1:X3 is identified:
             ("W:X1:X3" %in% decision  & 
                #all inconsistent interactions are one of the following:
                all( decision[is.na(match(decision, truth_INT))] %in% c("W:X1", "W:X3") )
             ) 
         )
      )
      |
        (metric=="interaction" & (truth_description == "3-way INT of W, X3, X2(4l)" |  #ordinal splitter 4 levels
                                    truth_description == "3-way INT of W, X3, X2(3l)" |  #ordinal splitter 3 levels
                                    truth_description == "3-way INT of W, X3, X2(l)") &  #cotinuous splitter
           #I.1.3.1. the same ignoring order
           ( setequal(Deci_vec, truth_INT) == T | 
               #I.1.3.2. as long as W:X1:X3 is identified:
               ("W:X2:X3" %in% decision  & 
                  #all inconsistent interactions are one of the following:
                  all( decision[is.na(match(decision, truth_INT))] %in% c("W:X2", "W:X3") )
               ) 
           )
        )
      |
        #----------------------------------------------------                                       
      #I.1.4. 2-way INT, two or three 2 way INTs
      #----------------------------------------------------
      (metric=="interaction" & (substr(truth_description,1,5) == "2-way" | #2-way INT
                                  truth_description == "Three 2-way INTs of W, X3, X1, X8" | #three 2-way INT
                                  truth_description == "Two 2-way INTs of W, X3, X1") #two 2-way INT
       & setequal(Deci_vec, truth_INT) ) 
    ), 
    1 
    , 
    ############################--------------------------       
    # I.2 "FALSE NEGATIVE" abuse term
    ############################--------------------------
    #SG/INT No. < truth     & at 1st glance, none of defined subgroups are the same as the truth;   & ALL SG/INT in the decision vector are nested in truth list, e.g. X3:X1 in X3:X1:X8   
    ifelse( N_deci < N_true  & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                        & all( unlist(DETECT_STRING(Truth_vec_L, Deci_vec)) ) == T, 2,  
            #SG/INT No. < truth     & at 1st glance, none of defined subgroups are the same as the truth;   & some or all of SG/INT in the decision vector are nested in truth list, e.g. X3:X1 in X3:X1:X8   
            ifelse( N_deci < N_true  & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                        & all( unlist(DETECT_STRING(Truth_vec_L, Deci_vec)) ) == F, 3,  
                    # subgroup No. < truth  & some of defined subgroups are the same as the true subgroups;         & some SG/INT in the truth list not included in the decision vector    
                    ifelse( N_deci < N_true  & length( Deci_vec[Deci_vec %in% Truth_vec]) > 0                        & (length( Truth_vec[Truth_vec  %in% Deci_vec]) < length(Truth_vec) ) , 4, 
                            ############################--------------------------       
                            # I.3 "FALSE POSITIVE" abuse term
                            ############################--------------------------       
                            # SG/INT No. > truth  & at 1st glance, none of SG/INT are the same as the truths;               & NONE of SG/INT in the truth vector are nested in decision list  (totally wrong)
                            ifelse( N_deci > N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == F, 5,  
                                    #SG/INT No. >= truth  & at 1st glance, none of SG/INT are the same as the truths;                & some or all of SG/INT in the truth vector (eg,"W:X1" "W:X3") are nested in decision list (eg, "W:X1:X10"    "W:X1:X10:X4" "W:X1:X10:X7"), e.g.
                                    ifelse( N_deci > N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == T, 6,  
                                            # SG/INT No. >= truth  & some SG/INT are the same as in the true subgroups;                      & some SG/INT in the decision list not included in the truth vector                         
                                            ifelse( N_deci > N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) > 0                        & ( length( Deci_vec[Deci_vec %in% Truth_vec]) < length(Deci_vec)  ) , 7, 
                                                    ############################--------------------------       
                                                    # I.4 
                                                    ############################-------------------------- 
                                                    # SG/INT No. = truth  & at 1st glance, none of SG/INT are the same as the truths;               & NONE of SG/INT in the truth vector are nested in decision list (totally wrong)
                                                    ifelse( N_deci == N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == F, 8,  
                                                            #SG/INT No. >= truth  & at 1st glance, none of SG/INT are the same as the truths;                & some or all of SG/INT in the truth vector (eg,"W:X1" "W:X3") are nested in decision list (eg, "W:X1:X10"    "W:X1:X10:X4"), e.g.
                                                            ifelse( N_deci == N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) == 0                       & all( unlist(DETECT_STRING(Deci_vec_L, Truth_vec)) ) == T, 9,  
                                                                    # SG/INT No. >= truth  & some SG/INT are the same as in the true subgroups;                      & some SG/INT in the decision list (eg "X3<=0 & X1<=0" "X3>0 & X1<=0" "X7<=-0.2 & X1>0" "X7>-0.2 & X1>0") not included in the truth vector (eg "X3<=0 & X1<=0" "X3<=0 & X1>0"  "X3>0 & X1<=0"  "X3>0 & X1>0" )                          
                                                                    ifelse( N_deci == N_true & length( Deci_vec[Deci_vec %in% Truth_vec]) > 0                        & ( length( Deci_vec[Deci_vec %in% Truth_vec]) < length(Deci_vec)  ) , 10, 
                                                                            11 )))))))))
    
  )
  
  #the following is for "subgroup" metric!!!
  sgpct_Nacc_Nalldeci <- length( Deci_vec[Deci_vec %in% Truth_vec])/length(Deci_vec)  #% identified true subgroups among identified subgroups
  sgpct_Nacc_Nalltrue <- length( Deci_vec[Deci_vec %in% Truth_vec])/length(Truth_vec) #% identified true subgroups among true subgroups, 
  
  return( list(value         = value,
               sg_Nacc_Nall  = sgpct_Nacc_Nalldeci , 
               sg_Nacc_NallT = sgpct_Nacc_Nalltrue))
}


#' OTHERTREE_DECI
#' 
#' Function that obtain subgroup decision from other tree-based method such as IT tree   
#' @param vote_subgroup the subgroup decision from voted tree structure
#' 
#' @return subgroup decision   
#' 
#' @export
OTHERTREE_DECI <- function(vote_subgroup){
  if(identical(vote_subgroup,"NA" ) ){
    Deci = "NA" 
  } else {
    Deci <- vote_subgroup$majority
  }
  return(Deci)
}



#' SPLITVAR_SIGN_SPLITVAL_LofDf_VT
#' 
#' Function that make subgroup defintion as a list of df, 
#' each df representing each condition of this subgroup definition and has 3 key variabs: parent_split_var_, parent_sign_, parent_split_val
#' @param Deci_VT_con0 subgroup definition string with "&" sign come from VirtualTWin 
#' 
#' @return subgroup defintion as a list of df, each df representing each condition of this subgroup definition   
#' 
#' @export

SPLITVAR_SIGN_SPLITVAL_LofDf_VT <- function( Deci_VT_con0 ) {
  #step0. remove "&" from original subgroup defintion to prepare
  Deci_VT_con <- stringr::str_split(Deci_VT_con0, " & ", simplify = T) #KEY
  
  #step1. make each subgroup defintion a list, each conditon is an element of the list
  Deci_VT_con_l <- as.list(Deci_VT_con)
  Deci_VT_con_replace1 <- lapply(Deci_VT_con_l,        function(x) stringr::str_replace_all(x, ">=", " >= "))
  Deci_VT_con_replace2 <- lapply(Deci_VT_con_replace1, function(x) stringr::str_replace_all(x, "< ", " < "))
  Deci_VT_con_split    <- lapply(Deci_VT_con_replace2, function(x) stringr::str_split(x, " ")) #KEY
  
  #step2. make a dataframe of this subgroup defintion to manipulate split variable, sign, and split value
  Deci_VT_con_split2 <- data.frame(matrix(unlist(Deci_VT_con_split), nrow=length(Deci_VT_con_split), byrow=TRUE)) %>%
    dplyr::rename(parent_split_var_ = X1, parent_sign_ = X2, parent_split_val_ = X3) %>%
    dplyr::mutate(parent_split_val_ = as.numeric(as.character(parent_split_val_) ) )
  #convert dataframe to list:
  Deci_VT_con_split3 <-  split(Deci_VT_con_split2, seq(nrow(Deci_VT_con_split2)))
  
  return(Deci_VT_con_split3)
}

#' SPLITVAR_SIGN_SPLITVAL_LofDf_CT 
#' 
#' Similar to SPLITVAR_SIGN_SPLITVAL_LofDf_VT, used in the VT2CF_format function below
#' Function that make subgroup defintion as a list of df, 
#' each df representing each condition of this subgroup definition and has 3 key variabs: parent_split_var_, parent_sign_, parent_split_val
#' @param Deci_CT_con0 subgroup definition string with "&" sign come from Causal foret or with it's format (VirtualTwin traformed into CF format)
#' 
#' @return subgroup defintion as a list of df, each df representing each condition of this subgroup definition   
#' 
#' @export

SPLITVAR_SIGN_SPLITVAL_LofDf_CT <- function( Deci_CT_con0 ) {
  #step0. remove "&" from original subgroup defintion to prepare
  Deci_CT_con <- stringr::str_split(Deci_CT_con0, " & ", simplify = T) #KEY
  
  #step1. make each subgroup defintion a list, each conditon is an element of the list
  Deci_CT_con_l <- as.list(Deci_CT_con)
  Deci_CT_con_replace1 <- lapply(Deci_CT_con_l,        function(x) stringr::str_replace_all(x, "<=", " <= "))
  Deci_CT_con_replace2 <- lapply(Deci_CT_con_replace1, function(x) stringr::str_replace_all(x, ">", " > "))
  Deci_CT_con_split    <- lapply(Deci_CT_con_replace2, function(x) stringr::str_split(x, " ")) #KEY
  
  #step2. make a dataframe of this subgroup defintion to manipulate split variable, sign, and split value
  Deci_CT_con_split2 <- data.frame(matrix(unlist(Deci_CT_con_split), nrow=length(Deci_CT_con_split), byrow=TRUE)) %>%
    dplyr::rename(parent_split_var_ = X1, parent_sign_ = X2, parent_split_val_ = X3) %>%
    dplyr::mutate(parent_split_val_ = as.numeric(as.character(parent_split_val_) ) )
  #convert dataframe to list:
  Deci_CT_con_split3 <-  split(Deci_CT_con_split2, seq(nrow(Deci_CT_con_split2)))
  
  return(Deci_CT_con_split3)
}


#' VT2CF_format
#' 
#' Function that transform the splitting format of VirtualTwin style to CausalForest style:
#' 1) binary/ordinal splitter >= n.5 to >n, < n.5 to <= n; 
#' 2) continuous splitter: >= n to > n-0.001, < n to <= n-0.001
#' @param Deci_VT_con0 subgroup definition string with "&" sign (may come from VirtualTWin or CF)
#' 
#' @return subgroup decision   
#' 
#' @export

VT2CF_format <- function(Deci_VT_con0){
  
  Deci_VT_con_split3 <- SPLITVAR_SIGN_SPLITVAL_LofDf_VT(Deci_VT_con0)
  #step3. adjust splitting format: 1) binary/ordinal splitter >= n.5 to >n, < n.5 to <= n; 2) continuous splitter: >= n to > n-0.001, < n to <= n-0.001
  Deci_VT_con_split4 <- lapply(Deci_VT_con_split3, 
                               function(x) 
                                 if( #scenario 1: binary/ordinal splitter >= 
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) %in% c(2:8) &
                                   x[1,"parent_sign_"]==">=" 
                                 ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = ">", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.5, split_val_round_posi ))
                                 } else if ( 
                                   #scenario 2: binary/oridinal splitter <   
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) %in% c(2:8) &
                                   x[1,"parent_sign_"]     =="<" 
                                 ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = "<=", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.5, split_val_round_posi ))
                                 } else if ( 
                                   #scenario 3: continuous splitter >=
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) > 8 &
                                   x[1,"parent_sign_"]==">=" ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = ">", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.001, split_val_round_posi ))
                                 } else if ( 
                                   #scenario 3: continuous splitter <
                                   length(unique (  eval(parse(text=paste("X$", x[1,"parent_split_var_"], sep = ""))) )) > 8 &
                                   x[1,"parent_sign_"]=="<" ) {
                                   x <- x %>% dplyr::mutate(parent_sign_ = "<=", 
                                                            parent_split_val_ = round(parent_split_val_ - 0.001, split_val_round_posi ))
                                 }
  )
  
  
  #step4. combine information from split variable, sign, and split value into one column 
  Deci_VT_con_split5 <- lapply(Deci_VT_con_split4, function(x) x %>% dplyr::mutate(condition = paste0(parent_split_var_, parent_sign_, parent_split_val_)) %>%
                                 dplyr::select(condition) %>%
                                 `colnames<-`( NULL )
  )
  
  #step5. collapse list of strings into one single string  
  Deci_VT_con_split6 <-  paste ( as.vector( unlist(Deci_VT_con_split5 ) ), collapse = " & " )
  
  return(Deci_VT_con_split6)
  
}



#' SG_CONDITION_2_INT
#' 
#' Function that convert each subgroup definition into interaction (e.g. a subgroup defined as X1<=0 & X3<=0 will be converted into W:X1:W3)
#' @param condition a specific subgroup definition (one group), e.g.  X1<=0 & X3<=0 
#' @return interaction, e.g. W:X1:W3  
#' 
#' @export

SG_CONDITION_2_INT <- function(condition){
  #make subgroup defintion as a list of df, each df representing each condition of this subgroup definition and has 3 key variabs: parent_split_var_, parent_sign_, parent_split_val
  splitvar_sign_splitval_L <-  SPLITVAR_SIGN_SPLITVAL_LofDf_CT(condition) 
  #combine conditions from a list of a subgroup into a dataframe 
  condition_key_var_df <- dplyr::bind_rows(splitvar_sign_splitval_L, .id = "column_label") 
  #extract key covariates in a vector
  condition_key_var_vec <- unique( str_sort( as.vector( dplyr::pull(condition_key_var_df, parent_split_var_)   ) ) )
  #make it into W:X1:X2 format
  condition_INT <- paste("W", paste(condition_key_var_vec, collapse = ":"), sep = ":")
  
  return(condition_INT)
  
}

#' SG2INT
#' 
#' Function that convert overall subgroup decision into interaction expressions 
#' e.g. subgroup decision as 1: X1<=0; 2: X1>0 will be converted into W:X1 interaction expression
#' @param condition a specific subgroup definition (one group), e.g.  X1<=0 & X3<=0 
#' @return interaction, e.g. W:X1:W3  
#' 
#' @export
SG2INT <- function(deci){
  if (identical(deci,"NA" )  ) {
    SG2INT_vec = "NA"
  } else if (identical(deci,"Unknown") ) {
    SG2INT_vec = "Unknown"
  } else {
    #extract the subgroup defintion column (e.g. X1<=0, X1>0)
    Deci_col <- deci %>% dplyr::select(subgroup) %>% as.data.frame() %>%  `colnames<-`( NULL ) 
    #make a list containing defintion of each subgroup as each elements of the list
    Deci_col_L <-  split(Deci_col, seq(nrow(Deci_col)))
    #convert each subgroup definition into interaction (e.g. a subgroup defined as X1<=0 & X3<=0 will be converted into W:X1:W3)
    Deci_SG_INT <- lapply(Deci_col_L, SG_CONDITION_2_INT )
    #combine all INT terms into one string, remove duplicate
    SG2INT_str_uniq <- paste(unique(unlist(Deci_SG_INT)), collapse = ' ')
    #turn string into vector
    SG2INT_vec <- strsplit(SG2INT_str_uniq, "\\s+")[[1]]
  }
  return(SG2INT_vec)
}




#' MSE_DltY
#' 
#' Function that extract MSE from MSE dataframe statify population and saved stratified population (subgroup) in a list. 
#' @param DltY_true the ture DltY 
#' @param DltY the dataframe of MSE obtained from function DltY_DATA
#' 
#' @return the list of MSE, residual by IPTW (ate) or SMR (att) weighting
#' 
#' @export
#' 

MSE_DltY <- function(DltY_true, DltY) {
  env=global_env()
  MSE_df <- dplyr::left_join(DltY_true, dplyr::select(DltY,  c(ID, SubgroupID, Definition,     ps,            iptw_u,     iptw_s,        smrw,
                                                               dltY_ate, dltY_ate_low, dltY_ate_up, dltY_att, dltY_att_low, dltY_att_up, dltY_crude, dltY_crude_low, dltY_crude_up)), 
                             by = "ID" ) %>%
    dplyr::filter(SubgroupID.y != "NA") %>% #when left_join with decision by aVirtualTwin, the total sample size < N (some population not included in the final decision, which leads to )
    dplyr::mutate(residual_ate =  dltY_ate_t -  dltY_ate,
                  mse_ate = mean(residual_ate^2),
                  residual_att =  dltY_att_t -  dltY_att,
                  mse_att = mean(residual_att^2))
  
  return(list(residual_ate  = MSE_df[1, "residual_ate"] , 
              MSE_ate       = MSE_df[1, "mse_ate"], 
              residual_att  = MSE_df[1, "residual_att"] , 
              MSE_att       = MSE_df[1, "mse_att"])
  )
}



#' CATE_MAX_SG
#' 
#' Function that identify subgroup with the largest CATE defined by: 
#' 1) binary outcome: Q(A)=(P(Y=1|T=1, X∈A) - P(Y=1|T=0, X∈A)) – (P(Y=1|T=1) - P(Y=1|T=0))
#' 2) continuous:     Q(A)=(E(Y=1|T=1, X∈A) - E(Y=1|T=0, X∈A)) – (E(Y=1|T=1) - E(Y=1|T=0))
#' @param DltY.data calculated dltY
#' 
#' @return the subgroup leading to the largest Q(A)
#' @export
#' 
CATE_MAX_SG <- function(DltY.data){
  
  CATE_L_true_te  <-   DltY.data %>%
    dplyr::select(SubgroupID, Definition, dltY_ate) %>%
    group_by(Definition) %>%
    summarise(dltY_ate = mean(dltY_ate), n=n()) %>%
    dplyr::mutate(Dlt_ATEvsCATE = (dltY_ate - ATE_all_te) ) %>%
    dplyr::mutate(Dlt_ATEvsCATE_abs = abs(dltY_ate - ATE_all_te) ) %>%
    dplyr::arrange(desc(Dlt_ATEvsCATE_abs))
  
  CATE_largestSG_true_te <- CATE_L_true_te$Definition[1]
  
  if ( CATE_largestSG_true_te=="W <100") {CATE_largestSG_true_te = "NA"}
  
  return(CATE_largestSG_true_te)
  
}

