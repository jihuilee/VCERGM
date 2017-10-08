#' Create permuted samples to calculate p-value
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param teststat Test statistic calculated based on observed networks
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param interior.knot Number of interior knots to create splines. Default is 10.
#' @param lambda.range Range of lambda (Tuning parameter). Default is seq(-3, 3, by = 0.1).
#' @param NPerm Number of permuted samples. Default is 100.
#' @param seed Seed number used to simulate bootstrap samples. Default is 12345.
#'
#' @importFrom splines bs
#' @importFrom ergm simulate.ergm
#' @export

permutation_test = function(object, networks, attr = NULL, teststat, Delta = NULL,
                            directed = FALSE, degree.spline = 3, interior.knot = 3,
                            lambda.range = seq(-3, 3, by = 0.1), NPerm = 100, seed = 12345)
{
  K = length(networks)

  Perm.seq = matrix(NA, nrow = NPerm, ncol = K)
  perm.teststat = rep(NA, NPerm)
  
  for(p in 1:NPerm)
  {
    set.seed(seed + p)
    # Permute the sequences to rearrange the observed networks
    Perm.seq[p, ] = perm.seq = sample(K)
    net.p = networks[perm.seq]
    Delta.p = Delta[perm.seq]

    # H1 (vcergm)
    vcergm1 = estimate_vcergm(object = object, networks = net.p, degree.spline = degree.spline, attr = attr,
                              lambda.range = lambda.range, interior.knot = interior.knot, directed = directed, constant = FALSE)

    # H0 (Constant)
    vcergm0 = estimate_vcergm(object, networks = net.p, degree.spline = degree.spline, attr = attr,
                              lambda.range = lambda.range, interior.knot = interior.knot, directed = directed, constant = TRUE)

    # Calculating test statistic
    cat("Calculating test statistic for permuted sample", p, "/", NPerm, "\n")
    perm.teststat[p] = test_statistic(object = object, networks = net.p, attr = attr,
                                      phi0 = vcergm0$phi.hat, phi1 = vcergm1$phi.hat, directed = directed, Delta = Delta.p)$teststat
  }
  pvalue =  sum(teststat < perm.teststat) / NPerm
  return(list(perm.teststat = perm.teststat, perm.pvalue = pvalue))

}
