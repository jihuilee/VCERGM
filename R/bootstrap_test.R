#' Create bootstrap samples to calculate p-value
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param phicoef0 Estimates of basis coefficients under H0
#' @param phi0 Estimates of phi(t) under H0
#' @param teststat Test statistic calculated based on observed networks
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param interior.knot Number of interior knots to create splines. Default is 10.
#' @param lambda.range Range of lambda (Tuning parameter). Default is seq(-3, 3, by = 0.1).
#' @param MCMC.burnin MCMC burnin sample size. Default is 10000.
#' @param MCMC.interval Interval between selected networks. Default is 1000. The first simulated network is the (MCMC.burnin + MCMC.interval). Thereafter, every (MCMC.interval)th network will be sample.
#' @param NBoot Number of bootstrap samples. Default is 1000.
#' @param seed Seed number used to simulate bootstrap samples. Default is 123.
#'
#' @importFrom splines bs
#' @importFrom ergm simulate.ergm
#' @export

bootstrap_test = function(object, networks, attr = NULL, phicoef0 = NULL, phi0 = NULL, teststat,
                          directed = FALSE, degree.spline = 3, interior.knot = 3,
                          lambda.range = seq(-3, 3, by = 0.1),
                          MCMC.burnin = 10000, MCMC.interval = 1000, NBoot = 1000, seed = 123)
{
  set.seed(seed)
  K = length(networks)
  Bnet = bs(1:K, df = degree.spline + 1 + interior.knot, degree = degree.spline, intercept = TRUE)

  #replacing object with current network formula
  z = deparse(object[[3]])

  if (is.null(phi0)) {phi0 = phicoef0 %*% t(Bnet)}

  # simulate bootstrap sample
  boot.networks = vector("list", NBoot)

  for (s in 1:K)
  {
    num.nodes = nrow(networks[[s]])
    Bs = matrix(Bnet[s,])

    if(is.vector(phi0) == TRUE) {coefs = as.vector(phi0[s])} else {coefs = as.vector(phi0[,s])}

    if (is.null(attr) == FALSE)
    {
      if (is.matrix(attr[[s]]) == FALSE) {
        attrs = vector("list", 1); attrs[[1]] = attr[[s]]; names(attrs) = "attr1"
      } else {attrs = vector("list", ncol(attr[[s]]))
      for (l in 1:ncol_attr) {attrs[[l]] = attr[[s]][,l]}
      names(attrs) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      nets = network(num.nodes, vertex.attr = attrs, directed = directed)
    } else{nets = network(num.nodes, directed = directed)}

    formula.s = as.formula(paste("nets ~ ", z, sep = ""))

    # Use an existing function in package 'ergm'
    sims = simulate(object = formula.s, coef = coefs, nsim = NBoot, seed = seed,
                    control = control.simulate(MCMC.burnin = MCMC.burnin, MCMC.interval = MCMC.interval))

    netarray = array(NA, dim = c(num.nodes, num.nodes, NBoot))

    if (NBoot == 1) { boot.networks[[1]][[s]] = as.matrix.network(sims, matrix.type = "adjacency") }

    if (NBoot > 1)
    { for (i in 1:NBoot) { boot.networks[[i]][[s]] = as.matrix(sims[[i]], matrix.type = "adjacency") } }
  }

  boot.teststat = rep(NA, NBoot)
  for (b in 1:NBoot)
  {

    net.b = boot.networks[[b]]

    suppressWarnings("glm.fit: fitted probabilities numerically 0 or 1 occurred")

    # H1
    vcergm.false = estimate_vcergm(object = object, networks = net.b, attr = attr,
                                   degree.spline = degree.spline, interior.knot = interior.knot,
                                   directed = directed,
                                   lambda.range = lambda.range, constant = FALSE)

    # H0 (Constant)
    vcergm.true = estimate_vcergm(object = object, networks = net.b, attr = attr,
                                  degree.spline = degree.spline, interior.knot = interior.knot,
                                  directed = directed,
                                  lambda.range = lambda.range, constant = TRUE)


    phi1 = as.matrix(vcergm.false$phi.hat)
    phi0 = as.matrix(vcergm.true$phi.hat)

    boot.teststat[b] = test_statistic(object = object, networks = net, attr = attr,
                                      phi0 = phi0, teststat = teststat[k],
                                      degree.spline = degree.spline, interior.knot = interior.knot,
                                      directed = directed, NBoot = NBoot)$boot.teststat

    cat("Calculating test statistic for bootstrap sample", b, "/", NBoot, "\n")
  }

  pvalue = sum(teststat < boot.teststat) / NBoot
  return(list(boot.networks = boot.networks, boot.teststat = boot.teststat, pvalue = pvalue))
}
