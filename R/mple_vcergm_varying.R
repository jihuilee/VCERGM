#' Calculate the MPLE of a specified VCERGM
#'
#' @param object A formula object of the form (network object) ~ <model terms>. Model terms take the same form of ERGM terms from R package 'ergm'.
#' @param networks A list of observed networks. It should have a list() object.
#' @param attr attr A list of vertex attributes. Default is NULL. (i.e. No attributes)
#' @param directed TRUE for analyzing directed networks. FALSE for analyzing undirected networks.
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

mple2 = function(object, networks, attr, directed, B, degree.spline, lambda.range, constant, Tol)
  {
  K = length(networks)
  missing.indx = which(sapply(networks, is.null))
  available.indx = setdiff(1:K, missing.indx)

  mean.num.nodes = floor(mean(unlist(lapply(networks, nrow)))) # Varying newtork size
  
  mean.denom = calculate_weight(object = object, nnodes = mean.num.nodes, attr = attr[[1]])
  
  # Save the result as list:
  # for each network in the time series, we calculate the change matrix and the vector of edges

  design = yy = ww = rep(list(NULL), length(available.indx))

  start1 = Sys.time()
  for (s in 1:length(available.indx)) {

    cat("Analyzing Network ", available.indx[s], "of ", K, "\n")

    #functions for time s
    Bu = matrix(B[s, ])
    #network at time point s
    nets = networks[[available.indx[s]]]

    #number of nodes at time s
    nnodes.s = nrow(nets)

    if (is.null(attr) == FALSE)
    {
      if (is.matrix(attr[[s]]) == FALSE) {
        attr.s = vector("list", 1); attr.s[[1]] = attr[[s]]; names(attr.s) = "attr1"
      } else {attr.s = vector("list", ncol(attr[[s]]))
      for (l in 1:ncol(attr[[s]])) {attr.s[[l]] = attr[[s]][,l]}
      names(attr.s) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      nets = network(nets, vertex.attr = attr.s, directed = directed)
    } else{nets = network(nets, directed = directed)}

    # replace object with current network formula
    z = deparse(object[[3]])
    formula.s = nonsimp_update.formula(object, nets ~ ., from.new = TRUE)

    # calculate the edges and the associated change matrix
    temp = ergmMPLE(formula.s, output = "matrix")

    # full network
    denom.s = calculate_weight(object = object, nnodes = nnodes.s, attr = attr[[s]])
    
    # observed edges for network
    yy[[s]] = temp$response

    # change matrix
    h.stats = temp$predictor / denom.s 

    # weights
    ww[[s]] = temp$weight
    #the (n choose 2) x pq part of design matrix

    design[[s]] = kronecker(t(Bu), h.stats)
  }
  end1 = Sys.time()

  # Unlist the elements and concatenate them.
  w = unlist(ww)
  y = unlist(yy)
  design.matrix = NULL
  for (i in 1:length(design)) {design.matrix = rbind(design.matrix, design[[i]])}

  # run the ergmMPLE once to get the coefficient names
  stat.names = unlist(strsplit(deparse(object[[3]]), " "))
#  stat.names = stat.names[!stat.names %in% c("+", "=", "TRUE)", "FALSE)")]
  stat.names = stat.names[!stat.names %in% c("+", "=", "", "TRUE", "T", "T)", "TRUE)", "FALSE", "F", "F)", "FALSE)", "diff")]

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

calculate_weight = function(object, nnodes, attr = NULL)
{
  fullnet = matrix(1, nrow = nnodes, ncol = nnodes)
  
  if (is.null(attr) == FALSE)
  {
    if (is.matrix(attr) == FALSE) {
      attr.s = vector("list", 1); attr.s = rep(1, nnodes); names(attr.s) = "attr1"
    } else {attr.s = vector("list", ncol(attr))
    for (l in 1:ncol(attr)) {attr.s[[l]] = rep(1, nnodes)}
    names(attr.s) = paste("attr", 1:ncol(attr), sep = "")
    }
    nets = network(fullnet, vertex.attr = attr.s, directed = directed)
  } else{nets = network(fullnet, directed = directed)}
  
  fullform = nonsimp_update.formula(object, nets ~ ., from.new = TRUE)
  return(c(ergmMPLE(fullform, output = "matrix")$predictor))
}
