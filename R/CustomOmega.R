#' Calculate Omega (Smoothness matrix)
#'
#' @param B A set of basis functions (K x q matrix)
#' @param available.indx A list of observed time points
#'
#' @importFrom splines bs

CustomOmega = function(B, available.indx) 
  {
  # Degree of spline
  d = attr(B, "degree")
  
  if (d < 2) {stop("The degree of Bspline should be at least 2. Please respecify.")}
  
  int = as.vector(attr(B, "knots")) # Interior knots
  boundary = attr(B, "Boundary.knots") # Boundary knots
  knots0 = sort(c(int, boundary))
  
  # Knots
  knots = c(rep(boundary[1], d + 1), int, rep(boundary[length(boundary)], d + 1))
  
  # Calculate the second derivatives
  B2 = bs(available.indx, knots = int, degree = d - 2, intercept = TRUE)

  # Matrix : T x q
  R = matrix(0, nrow = ncol(B), ncol = ncol(B) - 2)
  
  for (k in 3:nrow(R))
  {
    term1 = 1/(knots[k + d - 1] - knots[k])
    
    if (term1 == Inf) {terms = 1}
    R[(k - 2), (k - 2)] = term1 * 1/(knots[k + d - 1] - knots[k - 1])
    R[(k - 1), (k - 2)] = -term1 * (1/(knots[k + d - 1] - knots[k - 1]) + 1/(knots[k + d] - knots[k]))
    R[k, (k - 2)] = term1 * 1/(knots[k + d] - knots[k])
  }
  
  R = R * d * (d - 1)
  
  Omega0 =  R %*% t(B2)
  
  # Intervel lengths
  diff = diff(available.indx)
  
  # Integration Omega0 %*% t(Omega0) over (1,T)
  Omega = matrix(NA, nrow = nrow(Omega0), ncol = nrow(Omega0))
  for (i in 1:nrow(Omega0))
  {
    for (j in 1:nrow(Omega0))
    {
      vec = Omega0[i,] * Omega0[j,]
      sum = 0
      for (k in 1:(length(available.indx) - 1))
      {sum = sum + (vec[k] + vec[k + 1]) * diff[k] / 2}
      Omega[i, j] = sum
    }
  }
  
  return(Omega)
}