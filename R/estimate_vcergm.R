#' Estimate VCERGM coefficients
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'. 
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param interior.knot Number of interior knots to create splines. Default is 10.
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#' @param lambda.range Range of lambda (Tuning parameter). Default is seq(-3, 3, by = 0.1).
#' @param constant TRUE for constrained model of homogeneous VCERGM. FALSE for unconstrained model for heterogeneous VCERGM. Default is FALSE.
#' @param Tol Tolerance level used for calculating MPLE (IRLS iterations). Default is 0.01.
#'
#' @importFrom splines bs
#' @export

estimate_vcergm = function(object, networks, attr = NULL,
                           directed = c(TRUE, FALSE),
                           degree.spline = 3,
                           interior.knot = 10,
                           lambda.range = seq(-3, 3, by = 0.1), 
                           constant = FALSE, Tol = 0.01)
{
  
  # Check that input is OK
  # 1) Check that object is a formula object
  if (class(object) != "formula" & class(object) != "ergm") {
    stop("Argument object must be either a formula or ergm object. Please respecify.") }
  
  # 2) Check that networks is a list() object
  if (class(networks) != "list") {
    stop("Argument networks must be a list object. Please respecify.")
  }
  
  # Set basis and other required parameters
  # 1) Degrees of freedom of basis splines
  q = degree.spline + 1 + interior.knot 
  
  # 2) T: length of time series 
  K = length(networks)
  
  missing.indx = which(sapply(networks, is.null))
  available.indx = setdiff(1:K, missing.indx)
  
  if (q > K) {
    warning("The argument degrees of freedom is greater than the number of networks observed. 
            Respecifying this to be the number of networks observed - 1.")
    q = K - 1
  }
  
  #4) Adjust the input according to DEFAULT if not specified
  directed = directed[1]
  
  #5) Create basis function matrix B
  # this will be a T x q matrix
  B = bs(available.indx, df = q, degree = degree.spline, intercept = TRUE)
  Bfull = bs(1:K, df = q, degree = degree.spline, intercept = TRUE)
  
  #====================================================================
  #Estimate the coefficient matrix beta using MPLE
  
  # Varying network size
  num.nodes.K = unlist(lapply(networks, nrow))
  
#  if(length(unique(num.nodes.K)) == 1){
    reg = mple(object = object, networks = networks, attr = attr, directed = directed, B = B,
               degree.spline = degree.spline, lambda.range = lambda.range, constant = constant, Tol = Tol)
#  } else{
#    reg = mple2(object = object, networks = networks, attr = attr, directed = directed, B = B,
#               degree.spline = degree.spline, lambda.range = lambda.range, constant = constant, Tol = Tol)
#  }
  #ests = coef(reg$logistic.reg)
  ests = reg$phicoef
  phicoef.hat = matrix(ests, ncol = q)
  phi.hat = phicoef.hat %*% t(Bfull)
  lambda = reg$lambda
  rownames(phi.hat) = c(reg$stat.names)
  rownames(phicoef.hat) = c(reg$stat.names)
  
  time.list = reg$time.list
  time.mple = reg$time.mple

  # Return phicoef.hat and phi.hat
  # phicoef.hat is the functional coefficient matrix
  # phi.hat is the value of the statistics at each time point
  
  # time.list is the computing time for creating a list of elements for mple
  # time.mple is the computing time for mple
  
  result = list(phicoef.hat = phicoef.hat, phi.hat = phi.hat, lambda = lambda, 
                 time.list = time.list, time.mple = time.mple)
}