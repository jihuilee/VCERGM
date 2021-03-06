% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_vcergm.R
\name{simulate_vcergm}
\alias{simulate_vcergm}
\title{Simulate new networks at each time point}
\usage{
simulate_vcergm(object, attr = NULL, num.nodes.K, phi = NULL,
  phicoef = NULL, B = NULL, nsim = 100, MCMC.burnin = 10000,
  MCMC.interval = 1000, seed = 123, directed = c(TRUE, FALSE))
}
\arguments{
\item{object}{A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.}

\item{attr}{Node-specific attributes. Default is NULL.}

\item{num.nodes.K}{Number of nodes at each time point (K x 1 vector)}

\item{phi}{TRUE phi(t) functions (number of network statistics x K vector)}

\item{phicoef}{TRUE basis coefficients (p x q matrix)}

\item{B}{A set of basis functions (K x q matrix)}

\item{nsim}{Number of simulated networks. Default is 100.}

\item{MCMC.burnin}{MCMC burnin sample size. Default is 10000.}

\item{MCMC.interval}{Interval between selected networks. Default is 1000. The first simulated network is the (MCMC.burnin + MCMC.interval). Thereafter, every (MCMC.interval)th network will be sample.}

\item{seed}{Seed number used to simulate networks. Default is 123.}

\item{directed}{TRUE for analyzing directed networks. FALSE for analyzing undirected networks.}
}
\description{
Simulate new networks at each time point
}
