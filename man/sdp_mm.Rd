% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sdp_mm}
\alias{sdp_mm}
\title{SDP-MM method (convexified modularity maximization)}
\usage{
sdp_mm(S, lambda = as.numeric(c()), maxiter = 100L, tolerance = 0.01,
  admm_penalty = 1, verbose = 0L)
}
\arguments{
\item{S}{Input symmetric matrix}

\item{lambda}{Vector of possible values of SDP tuning parameter}

\item{maxiter}{Max number of iteration}

\item{tolerance}{ADMM tolerance}

\item{admm_penalty}{ADMM penalty parameter}

\item{verbose}{Level of verbosity}
}
\value{
A list containing SDP results
}
\description{
This function computes a solution path of SDP-MM
under a sequence of values of tuning parameter lambda
(see Convexified Modularity Maximization for Degree-corrected Stochastic 
Block Models by Yudong Chen,  Xiaodong Li, and Jiaming Xu for details)
}
