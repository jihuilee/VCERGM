#' Penalized logistic regression (GLM with roughness penalty on coefficient vector)
#'
#' @param y Response matrix of the penalized logistic regression
#' @param H Design matrix of the penalized logistic regression
#' @param weights Weights of each row of design matrix H
#' @param B A set of basis functions (K x q matrix)
#' @param available.indx A list of observed time points
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param constant TRUE for constrained model of homogeneous VCERGM. FALSE for unconstrained model for heterogeneous VCERGM. Default is FALSE.
#' @param lambda.range Range of lambda (Tuning parameter)
#' @param Tol Tolerance level used for calculating MPLE (IRLS iterations)
#' @export

penlogistic = function(y, H, weights, B, available.indx, degree.spline, constant, lambda.range, Tol)
  {
  # set initial values
  Niter = 100
  temp = dim(H)
  n = temp[1]
  q = ncol(B)
  p.q = temp[2]

  # set Omega (for evenly spaced time points)
  Omega = CustomOmega(B, available.indx)  
  
  Omega = kronecker(Omega, diag(1,p.q/q))
  
  niter = 0
  diff = 10000
  #initialize with logistic regression estimates
  
  # VCERGM
  if (constant == FALSE) {  
    logistic.reg = glm(y ~ H - 1, weights = weights, family = binomial) 
    phicoef = c(coef(logistic.reg))
  }
  
  # Flat function
  if (constant == TRUE) {
    H2 = H %*% kronecker(as.matrix(rep(1, q)), diag(1, p.q/q, p.q/q))
    logistic.reg2 = glm(y ~ H2 - 1, weights = weights, family = binomial) 
    phicoef0 = as.matrix(coef(logistic.reg2))
    phicoef = c(phicoef0 %*% t(as.matrix(rep(1, q))))
  }
  
  if (length(lambda.range) == 1) {
    if (lambda.range == 0) {
      return(list(phicoef = phicoef, lambda = 0))
    }
  }

    while (niter < Niter & diff > Tol) {
    phicoef_old = phicoef;
    print(diff)
    # critical values
    theta = H %*% phicoef
    eta = H %*% phicoef
    mu = db_fun(theta)
    
    # weight
    diag_W = weights/(ddb_fun(theta)*(dg_fun(mu))^2) # n x 1 vector of diagonal of weight W
    sw = sqrt(diag_W)
    
    # working response
    y_temp = H %*% phicoef + (y - mu) * dg_fun(mu)     
    y_working = sw * y_temp
    
    if (constant == FALSE) 
    {
      # working covariate
#      H_working = diag(c(sw)) %*% H
      H_working = vec.mat.product(c(sw), H)
      
      # tuning
      lambda = GCV(y_working, H_working,Omega, lambda.range)
      
      # update IRLS+penalty
      phicoef = solve(crossprod(H_working, H_working) + n*lambda*Omega) %*% crossprod(H_working, y_working)
    }
    
    if (constant == TRUE) 
    {
      # working covariate
#      H_working = diag(c(sw)) %*% H
      H_working = vec.mat.product(c(sw), H)
      term = kronecker(as.matrix(rep(1, q)), diag(1, p.q/q, p.q/q))
      
      # tuning
      lambda = GCV(y_working, H_working, Omega, lambda.range)
      
      H_working = diag(c(sw)) %*% H %*% term
      
      # update IRLS+penalty
      phicoef0 = solve(crossprod(H_working, H_working) + n*lambda * t(term) %*% Omega %*% term) %*% crossprod(H_working, y_working)
      phicoef = matrix(phicoef0 %*% t(as.matrix(rep(1, q))))
    }
    
    # stopping rule
    niter = niter + 1
    diff = norm(phicoef - phicoef_old, 'f') # Frobenius norm
  }
  return(list(phicoef = phicoef, lambda = lambda))
}

# entrywise functions for bernoulli glm
b_fun = function(theta){
  out = log(1 + exp(theta))
  return(out)
}

g_fun = function(mu){
  out = log(mu / (1 - mu))
  return(out)
}

ginv_fun = function(eta){
  out = exp(eta)/(1 + exp(eta))
  return(out)
}

db_fun = function(theta){
  out = exp(theta)/(1 + exp(theta))
  return(out)
}

ddb_fun = function(theta){
  out = exp(theta)/(1 + exp(theta))^2
  return(out)
}

dg_fun = function(mu){
  out = 1/(mu*(1 - mu))
  return(out)
}


# GCV
GCV = function(y, H, Omega, lambda.range){
  # this function select the best lambda for
  # ||y-H*phicoef||^2+n*lambda*phicoef'*Omega*phicoef 
  # through GCV (without refitting)
  # Note: y X Omega should be matrix
  
  temp = dim(H)
  n = temp[1]
  p = temp[2]
  
  # default lambda range
  lamrange = 10^(lambda.range) / n
  
  # record GCV score for each candidate lambda
  CV_score = rep(0, length(lamrange)) 
  
  for (i in 1:length(lamrange)) {
    lambda = lamrange[i]
    
    ############### H: length(y) x (q x p)
    ############### Omega: K x K
    
    yy = crossprod(y, y)
    HH = crossprod(H, H)
    yH = crossprod(y, H)
    
    InvMat = solve(HH + n*lambda*Omega) # p*p matrix   
    
    numerator = (1/n)*(yy - 2* yH %*% InvMat %*% t(yH) + t(y) %*% H %*% InvMat %*% HH %*% InvMat %*% t(yH))
    denominator = (1 - (1/n)*sum(diag(HH %*% InvMat)))^2
    CV_score[i] = numerator/denominator
  }
  ind = which.min(CV_score)
  bestlambda = lamrange[ind]
  return(bestlambda)
}

vec.mat.product = function(vec, mat)
{
  new = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  for(i in 1:length(vec))
  {
    new[i,] = vec[i] * mat[i,]
  }
  return(new)
}