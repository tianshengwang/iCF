% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iCF_SUBGROUP_DECISION.R
\name{DECI_treeN1_mCF_TF}
\alias{DECI_treeN1_mCF_TF}
\title{Function that check if root node of tree structure is accurate}
\usage{
DECI_treeN1_mCF_TF(
  decision,
  vote_D4_sg,
  vote_D3_sg,
  vote_D2_sg,
  vote_D4_t,
  vote_D3_t,
  vote_D2_t,
  tree_true
)
}
\arguments{
\item{decision}{the final decision obtained from the decision path}

\item{vote_D4_sg}{the list of D4 vote subgroup results}

\item{vote_D3_sg}{the list of D3 vote subgroup results}

\item{vote_D2_sg}{the list of D2 vote subgroup results}

\item{vote_D4_t}{the list of D4 vote tree results}

\item{vote_D3_t}{the list of D3 vote tree results}

\item{vote_D2_t}{the list of D2 vote tree results}

\item{tree_true}{the truth of tree structure}
}
\value{
if root node of tree structure is accurate
}
\description{
Function that check if root node of tree structure is accurate
}
