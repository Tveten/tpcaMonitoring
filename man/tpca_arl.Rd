% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arl_functions.R
\name{tpca_arl}
\alias{tpca_arl}
\title{Estimates average run length (ARL) for monitoring by tpca}
\usage{
tpca_arl(threshold, mu_x, Sigma_x, axes, n, w, n_sim)
}
\arguments{
\item{threshold}{A numeric specifying the threshold value for when a change
is declared.}

\item{mu_x}{A mean vector estimated from training data, representing the null distribution.}

\item{Sigma_x}{A covariance matrix estimated from training data, representing the null distribution.}

\item{n}{The number of observations to monitor for an estimate of the ARL.
See details.}

\item{w}{The window size. Number of recent time-points to consider for a
change.}

\item{n_sim}{The number of simulations to base the estimate on.}
}
\value{
An estimate of the average run length.
}
\description{
Description
}
\details{
Details
}
