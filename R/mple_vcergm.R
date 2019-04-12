#' Calculate the MPLE of a specified VCERGM
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
#' @param lag Lag in temporal networks. Default is 0 (i.e. VCERGM(0))
#' @param B A set of basis functions (K x q matrix)
#' @param degree.spline Degree of splines. Default is 3 (cubic splines).
#' @param lambda.range Range of lambda (Tuning parameter)
#' @param constant TRUE for constrained model of homogeneous VCERGM. FALSE for unconstrained model for heterogeneous VCERGM. Default is FALSE.
#' @param Tol Tolerance level used for calculating MPLE (IRLS iterations)
#'
#' @importFrom network network
#' @importFrom ergm ergmMPLE
#' @importFrom statnet.common nonsimp_update.formula
#' @export

mple = function(object, networks, attr, directed, lag, B, degree.spline, lambda.range, constant, Tol)
{
  K = length(networks)
  missing.indx = which(sapply(networks, is.null))
  available.indx = setdiff(1:K, missing.indx)

  # Save the result as list: For each network in the time series, we calculate the change matrix and the vector of edges
  networks2 = design = xx = yy = ww = rep(list(NULL), length(available.indx))

  start1 = Sys.time()
  for (s in 1:length(available.indx)) {

    cat("Analyzing Network ", available.indx[s], "of ", K, "\n")

    #functions for time s
    Bu = matrix(B[s, ])
    #network at time point s
    nets = networks[[available.indx[s]]]

    if (is.null(attr) == FALSE)
    {
      if (is.vector(attr[[s]])) {
        attr.s = vector("list", 1); attr.s[[1]] = attr[[s]]; names(attr.s) = "attr1"
      } else {attr.s = vector("list", ncol(attr[[s]]))
      for (l in 1:ncol(attr[[s]])) {attr.s[[l]] = attr[[s]][,l]}
      names(attr.s) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      nets = network(nets, vertex.attr = attr.s, directed = directed)
    } else{nets = network(nets, directed = directed)}

    networks2[[s]] = nets

    # replace object with current network formula
    formula.s = nonsimp_update.formula(object, nets ~ ., from.new = TRUE)

    # calculate the edges and the associated change matrix
    temp = ergmMPLE(formula.s, output = "matrix")

    # observed edges for network
    yy[[s]] = temp$response

    # change matrix
    xx[[s]] = temp$predictor

    # weights
    ww[[s]] = temp$weight
    #the (n choose 2) x pq part of design matrix

    if(lag == 0) {design[[s]] = kronecker(t(Bu), xx[[s]])}
  }
  end1 = Sys.time()

  stat.names = colnames(xx[[1]])

  theta = 0.001
  if(lag > 0)
  {
    h.stats = calculate_netstat(networks = networks, object = object, attr = attr)
    n.stats = nrow(h.stats)
    for(l in 1:length(networks))
    {
      if(l <= lag) {design[[l]] = kronecker(t(Bu), xx[[l]])} else{
        Bu = matrix(B[l, ])
        hstats.l = matrix(rep(h.stats[,l-1], nrow(xx[[l]])), ncol = n.stats, byrow = T)
        design[[l]] = kronecker(t(Bu), xx[[l]] - theta * hstats.l)
        }
    }
  }

  # Unlist the elements and concatenate them
  w = unlist(ww)
  y = unlist(yy)
  design.matrix = do.call("rbind", design)

  # run a penalized logistic regression of y on design.matrix to get pq x 1 parameter estimates
  # currently not using an intercept term
  q = dim(design.matrix)[2] / length(stat.names)

  start2 = Sys.time()
  logistic.reg = penlogistic(y = y, H = design.matrix, weights = w, B = B,
                             available.indx = available.indx, degree.spline = degree.spline,
                             constant = constant, lambda.range = lambda.range, Tol = Tol)
  end2 = Sys.time()

  return(list(phicoef = logistic.reg$phicoef, stat.names = stat.names, lambda = logistic.reg$lambda,
              time.list = end1 - start1, time.mple = end2 - start2))
}
