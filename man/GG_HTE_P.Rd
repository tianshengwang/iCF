% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GG_toolbox.R
\name{GG_HTE_P}
\alias{GG_HTE_P}
\title{Function that produce P-value for miCF}
\usage{
GG_HTE_P(df_D4, df_D3, df_D2, type_i, truth_describe, iflegend)
}
\arguments{
\item{df_D4}{depth 4 CF}

\item{df_D3}{depth 3 CF}

\item{df_D2}{depth 2 CF}

\item{type_i}{from all raw CF ("all") or those CF with best trees with the same structure as those majority voted best trees ("majority")}

\item{truth_describe}{describe the truth}

\item{iflegend}{if add legend}
}
\value{
the dataframe of residual & MSE
}
\description{
Function that produce P-value for miCF
}
