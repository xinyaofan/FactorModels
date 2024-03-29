% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdata.R
\name{sim1fact2}
\alias{sim1fact2}
\title{simulate data from 1-factor copula model}
\usage{
sim1fact2(n, parobj1, qcond1, copname1, ivect)
}
\arguments{
\item{n}{sample size}

\item{parobj1}{parameter vector of dimension d}

\item{qcond1}{function for copula conditional inverse cdf}

\item{copname1}{copula family in the copula model}

\item{ivect}{flag that is T if qcond1 have vectorized forms}
}
\value{
yy N*D simulate data in the u-scale

vv simulated latent variables
}
\description{
simulate data from 1-factor model
}
