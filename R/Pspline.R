#' Penalized B-spline function
#'
#' @param y Response matrix of the penalized logistic regression
#' @param H Design matrix of the penalized logistic regression
#' @param lambda.range Range of lambda (Tuning parameter)
#' @param B A set of basis functions (K x q matrix)
#' @param available.indx A list of observed time points
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).

Pspline = function(y, H, lambda.range, B, available.indx, degree.spline) 
  {
  # timevec of same length with y, indicating the time points when y is observed
  # y should be a vector of parameters estimated by ERGM at different time points (K*1)
  # H should be a K*q matrix, with each column corresp to a basis function evaluated at time points s=1,..,K
  # the output phicoef is a q*1 coefficient vector for the q basis functions. 
  # H %*% phicoef will be the fitted smooth curve for this ERGM parameter 
  
  # set Omega (for evenly spaced time points)
  Omega = CustomOmega(B, available.indx) 
  
  # use GCV to choose the best tuning parameter
  lambda = GCV(y, H, Omega, lambda.range)
  
  n = nrow(H)
  
  # update IRLS+penalty
  phicoef = solve(t(H) %*% H + n*lambda*Omega) %*% t(H) %*% y
  
  return(list(phicoef = phicoef, lambda = lambda))
}

