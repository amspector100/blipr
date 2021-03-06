% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cand_groups.R
\name{susie_groups}
\alias{susie_groups}
\title{Creates candidate groups based on a SuSiE model}
\usage{
susie_groups(alphas, X, q, max_pep = 1, max_size = 25, prenarrow = TRUE)
}
\arguments{
\item{alphas}{An L x p matrix of alphas from a SuSiE object.}

\item{X}{the N x p design matrix. If not NULL, will also add
hierarchical groups based on a correlation cluster of X.}

\item{q}{The nominal level at which to control the error rate}

\item{max_pep}{The maximum posterior error probability
(PEP) allowed in a candidate group. Default is 1.}

\item{max_size}{maximum allowable size for each group.}

\item{prenarrow}{If true, prenarrows the candidate groups
as described in the paper. Defaults to TRUE.}
}
\description{
Creates candidate groups based on a SuSiE model
}
