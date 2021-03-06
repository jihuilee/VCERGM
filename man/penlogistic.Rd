% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penalized_logistic.R
\name{penlogistic}
\alias{penlogistic}
\title{Penalized logistic regression (GLM with roughness penalty on coefficient vector)}
\usage{
penlogistic(y, H, weights, B, available.indx, degree.spline, constant,
  lambda.range, Tol)
}
\arguments{
\item{y}{Response matrix of the penalized logistic regression}

\item{H}{Design matrix of the penalized logistic regression}

\item{weights}{Weights of each row of design matrix H}

\item{B}{A set of basis functions (K x q matrix)}

\item{available.indx}{A list of observed time points}

\item{degree.spline}{Degree of splines. Default is 3 (cubic splines).}

\item{constant}{TRUE for constrained model of homogeneous VCERGM. FALSE for unconstrained model for heterogeneous VCERGM. Default is FALSE.}

\item{lambda.range}{Range of lambda (Tuning parameter)}

\item{Tol}{Tolerance level used for calculating MPLE (IRLS iterations)}
}
\description{
Penalized logistic regression (GLM with roughness penalty on coefficient vector)
}
