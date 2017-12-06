#' Hypothesis test for heterogeneity
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param interior.knot Number of interior knots to create splines. Default is 10.
#' @param phicoef0 Estimates of basis coefficients under H0
#' @param phicoef1 Estimates of basis coefficients under H1
#' @param phi0 Estimates of phi(t) under H0
#' @param phi1 Estimates of phi(t) under H1
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#'
#' @importFrom splines bs
#' @importFrom network network
#' @importFrom ergm summary.statistics
#' @importFrom ergm ergm.update.formula
#' @export

test_statistic = function(object, networks, attr = NULL, degree.spline = 3, interior.knot = 3,
                          phicoef0 = NULL, phicoef1 = NULL, phi0 = NULL, phi1 = NULL, directed = c(TRUE, FALSE), Delta = NULL)
{
  directed = directed[1]
  K = length(networks)
  B = bs(1:K, df = degree.spline + 1 + interior.knot, degree = degree.spline, intercept = TRUE)

  sum = 0
  if (is.null(Delta) == TRUE) {Delta = vector("list", length(networks))}
  for (s in 1:length(networks))
  {
    nets = networks[[s]]

    attrs = NULL

    if (is.null(attr) == FALSE)
    {
      if (is.matrix(attr[[s]]) == FALSE) {
        attrs = vector("list", 1); attrs[[1]] = attr[[s]]; names(attrs) = "attr1"
      } else {attrs = vector("list", ncol(attr[[s]]))
      for (l in 1:ncol(attr[[s]])) {attrs[[l]] = attr[[s]][,l]}
      names(attrs) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
    }

    if ((1 - is.null(phicoef0))*(1 - is.null(phicoef1)))
    {
      Bs = B[s,]
      stat.s = test.stat.s(nets, object, attrs, Bs = Bs,
                           phicoef0 = phicoef0, phicoef1 = phicoef1, directed = directed, Delta = Delta[[s]])
      temp = stat.s$teststat
      if (is.null(Delta[[s]]) == TRUE) {Delta[[s]] = stat.s$Delta}
    }

    if ((1 - is.null(phi0))*(1 - is.null(phi1)))
    {
      if (is.vector(phi0) == TRUE) {
        phi0s = phi0[s]
        phi1s = phi1[s]
      } else {
        phi0s = phi0[,s]
        phi1s = phi1[,s]
      }
      stat.s = test.stat.s(nets, object, attrs, phi0s = phi0s, phi1s = phi1s, directed = directed, Delta = Delta[[s]])
      temp = stat.s$teststat
      if (is.null(Delta[[s]]) == TRUE) {Delta[[s]] = stat.s$Delta}
    }

    sum = sum + temp
    cat("s is", s, "and sum is", sum, "\n")
  }
  return(list(teststat = as.numeric(sum), Delta = Delta))
}

test.stat.s = function(nets, object, attrs = NULL, Bs = NULL,
                       phicoef0 = NULL, phicoef1 = NULL, phi0s = NULL, phi1s = NULL, directed, Delta = NULL)
{
  z = deparse(object[[3]])
  netc = c(nets)

  # Network statistics
  stat = unlist(strsplit(deparse(object[[3]]), " "))
  stat = stat[!stat %in% c("+", "=", "TRUE)", "FALSE)")]
  nstat = length(stat)

  if (directed == FALSE) {
    netseq = netc[c(lower.tri(nets))]
  }

  if (directed == TRUE) {
    tri = as.logical(c(lower.tri(nets) + upper.tri(nets)))
    netseq = netc[tri]
  }

  if(is.null(Delta) == TRUE) {
    Delta = matrix(NA, nrow = length(netseq), ncol = nstat)

    for (i in 1:length(netseq))
    {
      net0 = net1 = netseq
      net0[i] = 0 # When X_ij = 0
      net1[i] = 1 # When X_ij = 1

      net00 = net11 = matrix(0, nrow = nrow(nets), ncol = ncol(nets))
      if (directed == FALSE) {
        net00[lower.tri(nets)] = net0
        net11[lower.tri(nets)] = net1

        # Make it symmetric
        net00 = net00 + t(net00)
        net11 = net11 + t(net11)
      }

      if (directed == TRUE) {
        tri = as.logical(lower.tri(nets) + upper.tri(nets))
        netseq = netc[tri]
        net00[tri] = net0
        net11[tri] = net1
      }

      if (is.null(attrs) == FALSE)
      {
        net00 = network(net00, vertex.attr = attrs, directed = directed)
        net11 = network(net11, vertex.attr = attrs, directed = directed)
      } else{
        net00 = network(net00, directed = directed)
        net11 = network(net11, directed = directed)
      }

      # When X_ij = 0
#      form0 = as.formula(paste("net00 ~ ", z, sep = ""))
      form0 = ergm.update.formula(object, net00 ~ .)

      # When X_ij = 1
#      form1 = as.formula(paste("net11 ~ ", z, sep = ""))
      form1 = ergm.update.formula(object, net11 ~ .)


      # Network statsitics when X_ij = 0 or 1
      g0 = summary.statistics(form0, directed = directed)
      g1 = summary.statistics(form1, directed = directed)

      Delta[i, ] = g1 - g0
    }
  }

if ((1 - is.null(phi0s))*(1 - is.null(phi1s)))
{
  term1 = (phi1s - phi0s) %*% t(Delta) %*% netseq
  term2 = sum(log((1 + exp(Delta %*% phi0s))/(1 + exp(Delta %*% phi1s))))
}

if ((1 - is.null(phicoef0))*(1 - is.null(phicoef1))*(1 - is.null(Bs)))
{
  term1 = t((phicoef1 - phicoef0) %*% Bs) %*% t(Delta) %*% netseq
  term2 = sum(log((1 + exp(Delta %*% (phicoef0 %*% Bs)))/(1 + exp(Delta %*% (phicoef1 %*% Bs)))))
}
  teststat = 2 * (term1 + term2)

  return(list(teststat = teststat, Delta = Delta))
}
