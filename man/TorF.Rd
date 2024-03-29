% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iCF_SUBGROUP_DECISION.R
\name{TorF}
\alias{TorF}
\title{function that check tree structure/subgroup are true or not}
\usage{
TorF(
  v_D4_tree,
  v_D4_tree_R,
  Desc_v_D4,
  v_D3_tree,
  v_D3_tree_R,
  Desc_v_D3,
  v_D2_tree,
  v_D2_tree_R,
  Desc_v_D2,
  truth.list
)
}
\arguments{
\item{v_D4_tree}{voted D4 tree}

\item{v_D4_tree_R}{voted D4 tree without leaves}

\item{Desc_v_D4}{decision from D4 tree}

\item{v_D3_tree}{voted D3 tree}

\item{v_D3_tree_R}{voted D3 tree without leaves}

\item{Desc_v_D3}{decision from D3 tree}

\item{v_D2_tree}{voted D2 tree}

\item{v_D2_tree_R}{voted D2 tree without leaves}

\item{Desc_v_D2}{decision from D2 tree}
}
\value{
subgroup in list format
}
\description{
function that check tree structure/subgroup are true or not
}
