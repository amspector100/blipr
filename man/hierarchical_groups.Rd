% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cand_groups.R
\name{hierarchical_groups}
\alias{hierarchical_groups}
\title{Creates hierarchically structured candidate groups
based on a distance matrix.}
\usage{
hierarchical_groups(
  samples,
  dist_matrix = NULL,
  X = NULL,
  max_pep = 1,
  max_size = 25,
  filter_sequential = FALSE
)
}
\arguments{
\item{samples}{(N,p)-shaped array of posterior samples
where a nonzero value indicates the presence of a signal.}

\item{dist_matrix}{A distance matrix corresponding to
distances between locations, used for hierarchical clustering.}

\item{X}{The design matrix in regression problems, which will
be used to create dist_matrix if dist_matrix is not provided.}

\item{max_pep}{The maximum posterior error probability
(PEP) allowed in a candidate group. Default is 1.}

\item{max_size}{maximum allowable size for each group.}

\item{filter_sequential}{If TRUE, ignore sequential groups
of variables to avoid duplication.}
}
\description{
Creates hierarchically structured candidate groups
based on a distance matrix.
}
