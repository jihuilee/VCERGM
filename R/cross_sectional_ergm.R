#' Calculate cross-sectional ERGM estimates
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param interior.knot Number of interior knots to create splines. Default is 10.
#' @param lambda.range Range of lambda (Tuning parameter). Default is seq(-3, 3, by = 0.1).
#'
#' @importFrom splines bs
#' @importFrom network network
#' @importFrom ergm ergm
#' @export

cross_sectional_ergm = function(object, networks, attr = NULL, directed = c(TRUE, FALSE),
                                degree.spline = 3, interior.knot = 10,
                                lambda.range = seq(-3, 3, by = 0.1))
  {
  directed = directed[1]

#  stat = unlist(strsplit(deparse(object[[3]]), " "))
#  stat = stat[!stat %in% c("+", "=", "TRUE)", "FALSE)")]
#  stat = stat[!stat %in% c("+", "=", "", "TRUE", "T", "T)", "TRUE)", "FALSE", "F", "F)", "FALSE)", "diff")]

#  hlength = length(stat)

  K = length(networks)
  missing.indx = which(sapply(networks, is.null))
  available.indx = setdiff(1:K, missing.indx)

  q = degree.spline + 1 + interior.knot
  B = bs(available.indx, df = q, degree = degree.spline, intercept = TRUE)

  Bfull = bs(1:K, df = q, degree = degree.spline, intercept = TRUE)

  # Cross-sectiona ERGMs
#  phi.hat = matrix(NA, hlength, K)
  phi.hat = NULL

  for (s in available.indx)
  {
    if (is.null(attr) == FALSE)
    {
      if (is.vector(attr[[s]])) {
        attr.s = vector("list", 1); attr.s[[1]] = attr[[s]]; names(attr.s) = "attr1"
      } else {attr.s = vector("list", ncol(attr[[s]]))
      for (l in 1:ncol(attr[[s]])) {attr.s[[l]] = attr[[s]][,l]}
      names(attr.s) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      sim.s = network(networks[[s]], vertex.attr = attr.s, directed = directed)
    } else {sim.s = network(networks[[s]], directed = directed)}

    z = deparse(object[[3]])
#    formula.s = as.formula(paste("sim.s ~ ", z, sep = ""))
#    formula.s = ergm.update.formula(object, sim.s ~ ., from.new = TRUE)
    formula.s = update(object, sim.s ~ .)

#    phi.hat[,s] = ergm(formula.s, estimate = "MPLE")$coef
    phi.hat = cbind(phi.hat, ergm(formula.s, estimate = "MPLE")$coef)
  }

  # Smoothing
#  phi.hat.smooth = matrix(NA, hlength, K)
  phi.hat.smooth = NULL

  for (k in 1:nrow(phi.hat))
  {
    smooth.coef = Pspline(y = phi.hat[k, available.indx], H = B, lambda.range = lambda.range,
                           available.indx = available.indx, B = B, degree.spline = degree.spline)$phicoef
#    phi.hat.smooth[k,] = Bfull %*% smooth.coef
    phi.hat.smooth = rbind(phi.hat.smooth, Bfull %*% smooth.coef)
  }

#  rownames(phi.hat) = rownames(phi.hat.smooth) = stat

  return(list(phi.hat = phi.hat, phi.hat.smooth = phi.hat.smooth))
}
