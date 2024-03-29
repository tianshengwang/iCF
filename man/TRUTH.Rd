% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_Truth_tree.R
\name{TRUTH}
\alias{TRUTH}
\title{This function MAKE TRUE TREE DATAFRAME FOR DIFFERENT SCENARIOS,
"tree_true" is less important "tree_true_subgroup" as different tree structures may lead to the same subgroup decision
We make the "tree_true" dataframe so that "tree_true_subgroup" can be generated
bo4= binary:ordinary with 4 level
bo3= binary:ordinal with 3 level
1o3= 1 ordinal 3 levels
1o4= 1 ordinal 4 levels}
\usage{
TRUTH(intTRUE)
}
\arguments{
\item{intTRUE}{scenario name which provide information on n-way/antigonism or synergistic/moderate or strong/binary, categorical, or continuous, e.g."2W_Ant_Mo_1b"}
}
\value{
The subsetted list that all are identical to reference list.
}
\description{
This function MAKE TRUE TREE DATAFRAME FOR DIFFERENT SCENARIOS,
"tree_true" is less important "tree_true_subgroup" as different tree structures may lead to the same subgroup decision
We make the "tree_true" dataframe so that "tree_true_subgroup" can be generated
bo4= binary:ordinary with 4 level
bo3= binary:ordinal with 3 level
1o3= 1 ordinal 3 levels
1o4= 1 ordinal 4 levels
}
