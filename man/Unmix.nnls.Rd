% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Unmixing_Methods.R
\name{Unmix.nnls}
\alias{Unmix.nnls}
\title{Non-negative least squares unmixing}
\usage{
Unmix.nnls(dat, M, bg = NULL, n_cores = 1)
}
\arguments{
\item{dat}{dataframe, data to be unmixed}

\item{M}{unmixing matrix}

\item{bg}{vector of constant background to be subtracted from each event}

\item{n_cores}{Number of cores to use during unmixing (default: 1)}
}
\value{
dataframe with component abundances
}
\description{
Non-negative least squares unmixing
}
